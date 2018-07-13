#!/usr/bin/env python

# (C) 2018 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

import tempfile

from cuberecord.ports import (
    RawRecordInputPort,
    RecordOutputPort,
    ShardOutputPort,
    ShardInputPort,
    ShardCollectionInputPort,
    ShardCollectionOutputPort,
)

from orionclient.session import APISession
from orionclient.types import ShardCollection, Shard

from datarecord import OEReadRecords
from cuberecord import TypedShard
from cuberecord.record_types import ShardTypes

from floe.api.cubes import SourceCube, ComputeCube, ParallelMixin
from floe.api.parameter import (
    StringParameter,
    IntegerParameter,
    BooleanParameter,
)


class RecordsToRecordCollectionConverter(ComputeCube):

    title = "Record to Record Collection"
    description = """
A writer that writes a stream of records to a Collection of record Shards (i.e. oedb).
                  """

    collection_setup = ShardCollectionInputPort("collection_setup", initializer=True)
    intake = RawRecordInputPort("intake")
    success = ShardOutputPort("success")

    def begin(self):
        for col in self.collection_setup:
            self.log.debug("Writing to Collection {} {}", col.id, col.name)
            self.collection = col
            break

    def process(self, record, port):
        with tempfile.NamedTemporaryFile(suffix=".oedb") as temp:
            ofs = open(temp.name, "wb")
            ofs.write(record)
            ofs.close()
            shard = TypedShard.create(
                collection=self.collection, shard_type=ShardTypes.Record
            )
            shard.upload_file(temp.name)
        self.log.debug("Added Shard {} to Collection {} {}", shard.id, self.collection.id, self.collection.name)
        self.success.emit(shard)


class RecordsToRecordCollectionConverterParallel(ParallelMixin, RecordsToRecordCollectionConverter):

    title = "Parallel " + RecordsToRecordCollectionConverter.title


class CollectionCompletionCube(ComputeCube):

    title = "Collection Completion"

    description = """
A sink cube that saves Shards and optionally close or delete the Collection.
    """

    collection_setup = ShardCollectionInputPort("collection_setup", initializer=True)
    intake = ShardInputPort("intake")

    save = BooleanParameter(
        "save", default=True, help_text="False will delete the Collection"
    )
    close = BooleanParameter(
        "close",
        default=False,
        help_text="Closing a collection means no shards can be added",
    )

    def begin(self):
        for col in self.collection_setup:
            self.log.info("Completing Collection {} {}", col.id, col.name)
            self.collection = col

    def process(self, shard, port):
        if self.args.save:
            shard.close()

    def end(self):
        if not self.args.save:
            APISession.delete_resource(self.collection)
        elif self.args.close:
            self.collection.close()


class CollectionCreationCube(SourceCube):

    title = "Collection Creation"

    description = """
An Collection initializer cube.
    """

    collection_setup = ShardCollectionOutputPort("collection_setup")

    collection_name = StringParameter("collection_name", required=True)

    def __iter__(self):
        self.collection = ShardCollection.create(
            APISession,
            name=self.args.collection_name
        )
        yield self.collection


class CollectionReader(SourceCube):

    title = "Collection Reader"
    description = "Reads a Collection and emits Shards"
    success = ShardOutputPort("success")

    collection = IntegerParameter(
        "collection", required=True, title="Collection to use as input"
    )
    shard_ids = StringParameter(
        "shard_ids",
        default="",
        title="Comma delimited ranges of shard ids (e.g. '1-3,5-7')",
    )
    limit = IntegerParameter(
        "limit",
        default=0,
        min_value=0,
        title="How many shards should be output, 0 indicates all within the library",
    )
    offset = IntegerParameter(
        "offset",
        default=0,
        min_value=0,
        title="How many shards should be skipped before emitting",
    )

    def begin(self):
        self.collection = APISession.get_resource(
            ShardCollection,
            self.args.collection
        )

    def __iter__(self):
        count = 0
        if not self.args.shard_ids:
            for shard in self.collection.list_shards():
                if self.args.offset and count < self.args.offset:
                    continue
                if self.args.limit and count >= self.args.limit:
                    break
                count += 1
                yield shard
        else:
            # https://stackoverflow.com/questions/6405208/how-to-convert-numeric-string-ranges-to-a-list-in-python
            # TODO: turn into a generator
            def parse_ranges(x):
                result = []
                for part in x.split(","):
                    if "-" in part:
                        a, b = part.split("-")
                        a, b = int(a), int(b)
                        result.extend(range(a, b + 1))
                    else:
                        a = int(part)
                        result.append(a)
                return result

            for shard_id in parse_ranges(self.args.shard_ids):
                if self.args.offset and count < self.args.offset:
                    continue
                if self.args.limit and count >= self.args.limit:
                    break
                count += 1
                yield Shard(collection=self.collection, id=shard_id)


class RecordShardToRecordConverter(ComputeCube):

    title = "Record Shard to Record Converter"
    description = "Reads a Shard of records and converts to records"
    intake = ShardInputPort("intake")
    success = RecordOutputPort("success")

    def process(self, shard, port):
        self.log.info("Downloading Shard {} in Collection {}", shard.id, shard.collection)
        with tempfile.NamedTemporaryFile(suffix=".oedb") as temp:
            shard.download_to_file(temp.name)
            with open(temp.name, "rb") as ifs:
                self.log.info("Reading Shard {} in Collection {}", shard.id, shard.collection)
                for record in OEReadRecords(ifs):
                    self.success.emit(record)


class RecordsShardToRecordConverterParallel(ParallelMixin, RecordShardToRecordConverter):

    title = "Parallel " + RecordShardToRecordConverter.title

from floe.api import WorkFloe
from floe.api.hyper import HyperCube

from oecubeutils.cubes import (
    CollectionReader, RecordsShardToRecordConverterParallel
)


wf = WorkFloe('Shard Collection Reader')

reader = CollectionReader('Collection Reader')
reader.promote_parameter(
    'collection',
    promoted_name='collection',
    description="Identifier of Collection, recommended to use ocli to find collections"
)
converter = RecordsShardToRecordConverterParallel('Converter')

wf.add_cubes(reader, converter)

reader.success.connect(converter.intake)

COllectionReader = HyperCube(wf)

COllectionReader.promote_port(converter, 'success')

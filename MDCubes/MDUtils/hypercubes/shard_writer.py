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

# from floe.api import WorkFloe
# from floe.api.hyper import HyperCube

# from MDCubes.MDUtils.cubes import (CollectionCreationCube,
#                                    CollectionCompletionCube,
#                                    RecordsToRecordCollectionConverterParallel)

# wf = WorkFloe('Shard Collection writer')

# creator = CollectionCreationCube('Collection Creator')
# creator.promote_parameter('collection_name', promoted_name='collection_name')
# writer = RecordsToRecordCollectionConverterParallel('writer')
# completer = CollectionCompletionCube('Collection Completer')

# wf.add_cubes(creator, completer, writer)

# creator.collection_setup.connect(writer.collection_setup)
# creator.collection_setup.connect(completer.collection_setup)
# writer.success.connect(completer.intake)

# CollectionWriter = HyperCube(wf)

# CollectionWriter.promote_port(writer, 'intake')

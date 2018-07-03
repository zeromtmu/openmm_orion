from floe.api import WorkFloe
from floe.api.hyper import HyperCube

from oecubeutils.cubes import (
    CollectionCreationCube, CollectionCompletionCube,
    RecordsToRecordCollectionConverterParallel
)


wf = WorkFloe('Shard Collection writer')

creator = CollectionCreationCube('Collection Creator')
creator.promote_parameter('collection_name', promoted_name='collection_name')
writer = RecordsToRecordCollectionConverterParallel('writer')
completer = CollectionCreationCube('Collection Completer')

wf.add_cubes(creator, completer, writer)

creator.collection_setup.connect(writer.collection_setup)
creator.collection_setup.connect(completer.collection_setup)
writer.success.connect(completer.intake)

HyperCycleCube = HyperCube(wf)

HyperCycleCube.promote_port(writer, 'intake')

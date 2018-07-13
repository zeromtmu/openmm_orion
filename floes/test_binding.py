from floe.api import WorkFloe

from YankCubes.cubes import YankBindingFECube

from cuberecord import (DataSetWriterCube,
                        DataSetReaderCube)

job = WorkFloe('Binding Affinity')

job.description = """
Testing
"""

job.classification = [['BindingFreeEnergy', 'Yank']]
job.tags = [tag for lists in job.classification for tag in lists]


isys = DataSetReaderCube("SystemReader", title="System Reader")
isys.promote_parameter("data_in", promoted_name="system", title='System Input File',
                       description="System file name")

abfe0 = YankBindingFECube("ABFE0", title="ABFE")

ofs = DataSetWriterCube('ofs', title='Out')
ofs.promote_parameter("data_out", promoted_name="out")

fail = DataSetWriterCube('fail', title='Failures')
fail.set_parameters(data_out='fail.oedb')


job.add_cubes(isys, abfe0, ofs, fail)
isys.success.connect(abfe0.intake)
abfe0.success.connect(ofs.intake)
abfe0.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()

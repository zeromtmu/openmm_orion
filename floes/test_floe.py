from floe.api import WorkFloe
from cuberecord import DataSetWriterCube
from LigPrepCubes.ports import LigandReaderCube
from ProtPrepCubes.ports import ProteinReaderCube
from LigPrepCubes.cubes import LigandChargeCube
from ComplexPrepCubes.cubes import ComplexPrepCube
from ComplexPrepCubes.cubes import HydrationCube
# from ComplexPrepCubes.cubes import SolvationCube
from ForceFieldCubes.cubes import ForceFieldCube
from MDCubes.OpenMMCubes.cubes import OpenMMminimizeCube
from MDCubes.OpenMMCubes.cubes import OpenMMNvtCube
from MDCubes.OpenMMCubes.cubes import OpenMMNptCube

# Declare and document floe
floe = WorkFloe('Testing Workfloe')
floe.description = "Testing Ligand Reader"
floe.classification = [["Ligand Preparation"]]
floe.tags = ["OpenEye"]

# Declare Cubes
ifs = LigandReaderCube('LigandReader')
lig_charge = LigandChargeCube("LigandCharge")
prot_reader = ProteinReaderCube("ProteinReader")
complx = ComplexPrepCube("ComplexPrep")
hydration = HydrationCube("Hydration")
# hydration = SolvationCube("Solvation")
ff = ForceFieldCube("ForceField")
minimize = OpenMMminimizeCube("Minimization")
nvt = OpenMMNvtCube("NVT")
npt = OpenMMNptCube("NPT")

ofs = DataSetWriterCube('ofs')

# Add cubes to floe
floe.add_cube(ifs)
floe.add_cube(lig_charge)
floe.add_cube(prot_reader)
floe.add_cube(complx)
floe.add_cube(hydration)
floe.add_cube(ff)
floe.add_cube(minimize)
floe.add_cube(nvt)
floe.add_cube(npt)
floe.add_cube(ofs)

# Promote parameters
ifs.promote_parameter('data_in',
                      promoted_name='ligands',
                      title='Ligand Reader')
prot_reader.promote_parameter('data_in',
                              promoted_name='protein',
                              title='Protein reader')

ofs.promote_parameter('data_out',
                      promoted_name='out',
                      title='Output File of Molecules')

# Cube Connections
ifs.success.connect(lig_charge.intake)
lig_charge.success.connect(complx.intake)
prot_reader.success.connect(complx.protein_port)
complx.success.connect(hydration.intake)
hydration.success.connect(ff.intake)
ff.success.connect(minimize.intake)
minimize.success.connect(nvt.intake)
nvt.success.connect(npt.intake)
npt.success.connect(ofs.intake)

if __name__ == "__main__":
    floe.run()


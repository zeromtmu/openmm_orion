from floe.api import (parameter,
                      ParallelMixin)

from cuberecord import (OERecordComputeCube,
                        OEField)
from cuberecord.constants import DEFAULT_MOL_NAME
from datarecord import (Types,
                        Meta,
                        ColumnMeta)

import traceback
from OpenMMCubes.utils import ParmedData
from oeommtools.utils import split
from TrjAnalysis.utils import BoundingBox
import sstmap as sm

from openeye import oechem
import mdtraj as md

class SSTMapSetCube(ParallelMixin, OERecordComputeCube):
    version = "0.0.0"
    title = "Water Thermodynamics by Using SSTMap"
    description = """
    SSTMap performs Water Thermodynamics analysis.
    SSTMaps supports hydration site analysis (HSA) 
    and Grid Inhomogeneous Solvation Theory (GIST). 
    
    SSTMap has been developed at Kurtzman Lab Lehman College 
    For more details, please visit 
    sstmap.org @ https://github.com/KurtzmanLab/SSTMap
    """
    classification = ["SSTMap Analysis"]
    tags = [tag for lists in classification for tag in lists]

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 43200},  # Default 12 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    center_ligand = parameter.BooleanParameter(
        'center_ligand',
        default=True,
        help_text='Center the grid around the ligand instead of the protein'

    )

    trj_fn = parameter.StringParameter(
        'trj_fn',
        default='trj.nc',
        help_text='Trajectory file name'
    )

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):

        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)

            field_system = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol,
                                   meta=ColumnMeta().set_option(Meta.Hints.Chem.PrimaryMol))

            if not record.has_value(field_system):
                self.log.warn("Missing molecule '{}' field".format(field_system.get_name()))
                self.failure.emit(record)
                return

            system = record.get_value(field_system)

            # Split the complex in components in order to apply the FF
            protein, ligand, water, excipients = split(system, ligand_res_name='LIG')

            self.log.info("\nProtein atom numbers = {}\nLigand atom numbers = {}\n"
                          "Water atom numbers = {}\nExcipients atom numbers = {}".format(protein.NumAtoms(),
                                                                                         ligand.NumAtoms(),
                                                                                         water.NumAtoms(),
                                                                                         excipients.NumAtoms()))

            if opt['center_ligand'] and ligand.NumAtoms():
                bb = BoundingBox(ligand, save_pdb=True, scale_factor=1.5)
            elif protein.NumAtoms():
                bb = BoundingBox(protein, save_pdb=True, scale_factor=1.0)
            else:
                oechem.OEThrow.Error("Protein and Ligand molecules have not been detected")

            field_parmed = OEField("Parmed", ParmedData)

            if not record.has_value(field_parmed):
                self.log.warn("Missing molecule '{}' field".format(field_parmed.get_name()))
                self.failure.emit(record)

            parmed_structure = record.get_value(field_parmed)

            system_top_fn = "system.prmtop"

            parmed_structure.save(system_top_fn, overwrite=True)
            parmed_structure.save("system.pdb", overwrite=True)

            center = (bb[0]+bb[1])*0.5

            trj = md.load(opt['trj_fn'], top=system_top_fn)

            gist = sm.GridWaterAnalysis(
                system_top_fn,
                opt['trj_fn'],
                start_frame=0, num_frames=trj.n_frames,
                grid_center=center,
                grid_dimensions=[bb[1][0]-bb[0][0], bb[1][1]-bb[0][1], bb[1][2]-bb[0][2]],
                prefix="phen")

            gist.print_system_summary()
            #gist.calculate_grid_quantities()

            # Emit the ligand
            self.success.emit(record)

        except Exception as e:
            # Attach an error message to the molecule that failed
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return

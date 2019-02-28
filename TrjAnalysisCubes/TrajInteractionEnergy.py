import traceback

from floe.api import ParallelMixin, parameter

from cuberecord import OERecordComputeCube

from datarecord import (OEField,
                        OERecord,
                        Types,
                        OEFieldMeta,
                        Meta)

from Standards import Fields

import TrjAnalysisCubes.utils as utl
import TrjAnalysisCubes.TrajMMPBSA_utils as mmpbsa


class TrajInteractionEnergyCube(ParallelMixin, OERecordComputeCube):
    title = "Trajectory Interaction Energies"
    version = "0.0.0"
    classification = [["Energy Analysis"]]
    tags = ['OEChem', 'OpenMM', 'TrajAnalysis', 'MMPBSA']
    description = """
    Protein-ligand interaction energies are calculated on an existing MD trajectory.
    The trajectory is taken from pre-existing protein and ligand trajectory OEMols.
    The forcefield used is taken from the parmed object associated with the trajectory
    OEMols.

    Input:
    -------
    Data Record with the ligand and Protein trajectory OEMols; the trajectory frames are
    included as conformers on the molecule. The associated parmed object is also on the
    record.

    Output:
    -------
    Data Record - The various energy components associated with the protein-ligand
    interaction energies are attached to the record as per-frame vectors of floats.
    This includes the MM interaction potential energies and their ligand, protein
    and complex components. The energy units are in kcal/mol.
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            opt = self.opt
            # Logger string
            opt['Logger'].info(' Beginning TrajInteractionEnergyCube')
            system_title = utl.RequestOEFieldType( record, Fields.title)
            opt['Logger'].info('{} Attempting to compute MD Traj protein-ligand Interaction energies'
                .format(system_title) )

            # Check that the OETraj analysis has been done
            analysesDone = utl.RequestOEField( record, 'AnalysesDone', Types.StringVec)
            if 'OETraj' not in analysesDone:
                raise ValueError('{} does not have OETraj analyses done'.format(system_title) )
            else:
                opt['Logger'].info('{} found OETraj analyses'.format(system_title) )

            # Extract the relevant traj OEMols from the OETraj record
            oetrajRecord = utl.RequestOEField( record, 'OETraj', Types.Record)
            opt['Logger'].info('{} found OETraj record'.format(system_title) )
            ligTraj = utl.RequestOEField( oetrajRecord, 'LigTraj', Types.Chem.Mol)
            opt['Logger'].info('{} #atoms, #confs in ligand traj OEMol: {}, {}'
                .format( system_title, ligTraj.NumAtoms(), ligTraj.NumConfs()) )
            protTraj = utl.RequestOEField( oetrajRecord, 'ProtTraj', Types.Chem.Mol)
            opt['Logger'].info('{} #atoms, #confs in protein traj OEMol: {}, {}'
                .format( system_title, protTraj.NumAtoms(), protTraj.NumConfs()) )

            # Extract the parmed object from the parent record
            prmed = utl.RequestOEFieldType( record, Fields.pmd_structure)

            # Compute interaction energies for the protein, ligand, and complex subsystems
            intE, cplxE, protE, ligE = mmpbsa.ProtLigInteractionEFromParmedOETraj(
                                       prmed, ligTraj, protTraj)
            if intE is None:
                raise ValueError('{} Calculation of Interaction Energies failed'.format(system_title) )

            # protein and ligand traj OEMols now have parmed charges on them; save these
            oetrajRecord.set_value(OEField('LigTraj', Types.Chem.Mol), ligTraj)
            oetrajRecord.set_value(OEField('ProtTraj', Types.Chem.Mol), protTraj)
            record.set_value(OEField('OETraj', Types.Record), oetrajRecord)

            # Create new record with traj interaction energy results
            opt['Logger'].info('{} writing trajIntE OERecord'.format(system_title) )
            trajIntE = OERecord()
            #
            intE_field = OEField("protein_ligand_interactionEnergy", Types.FloatVec,
                                meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajIntE.set_value( intE_field, intE)
            #
            ligE_field = OEField("ligand_intraEnergy", Types.FloatVec,
                                meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajIntE.set_value( ligE_field, ligE)
            #
            protE_field = OEField("protein_intraEnergy", Types.FloatVec,
                                 meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajIntE.set_value( protE_field, protE)
            #
            cplxE_field = OEField("complex_intraEnergy", Types.FloatVec,
                                 meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajIntE.set_value( cplxE_field, cplxE)
            # Add the trajIntE record to the parent record
            record.set_value( OEField( 'TrajIntE', Types.Record), trajIntE)
            analysesDone.append( 'TrajIntE')
            record.set_value( OEField( 'AnalysesDone', Types.StringVec), analysesDone)
            opt['Logger'].info('{} finished writing trajIntE OERecord'.format(system_title) )

            self.success.emit(record)

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            opt['Logger'].info('Exception {} in TrajInteractionEnergyCube on {}'.format(str(e),system_title) )
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


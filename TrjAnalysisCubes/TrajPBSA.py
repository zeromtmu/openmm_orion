import traceback

from floe.api import ParallelMixin

from cuberecord import OERecordComputeCube

from datarecord import (OEField,
                        OERecord,
                        Types,
                        OEFieldMeta,
                        Meta)

from Standards import Fields

import TrjAnalysisCubes.utils as utl

import TrjAnalysisCubes.TrajMMPBSA_utils as mmpbsa

import numpy as np

class TrajPBSACube(ParallelMixin, OERecordComputeCube):
    title = "Trajectory Poisson-Boltzmann and Surface Area Energies"
    version = "0.0.0"
    classification = [["Energy Analysis"]]
    tags = ['OEChem', 'Zap', 'TrajAnalysis', 'MMPBSA']
    description = """
    Protein-ligand interaction solvation energies are calculated on an existing MD trajectory.
    The trajectory is taken from pre-existing protein and ligand trajectory OEMols.
    The Poisson-Boltzmann and Surface Area methods in the OEZap toolkits are  used.

    Input:
    -------
    Data Record with the ligand and Protein trajectory OEMols; the trajectory frames are
    included as conformers on the molecule.

    Output:
    -------
    Data Record - The various energy components associated with the Poisson-Boltzmann and
    Surface Area energies are attached to the record as per-frame vectors of floats.
    The energy units are in kcal/mol.
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
            opt['Logger'].info(' Beginning TrajPBSACube')
            system_title = utl.RequestOEFieldType(record, Fields.title)
            opt['Logger'].info('{} Attempting to compute MD Traj PBSA energies'.format(system_title))

            # Check that the OETraj analysis has been done
            analysesDone = utl.RequestOEField( record, 'AnalysesDone', Types.StringVec)
            if 'OETraj' not in analysesDone:
                raise ValueError('{} does not have OETraj analyses done'.format(system_title))
            else:
                opt['Logger'].info('{} found OETraj analyses'.format(system_title))

            # Extract the relevant traj OEMols from the OETraj record
            oetrajRecord = utl.RequestOEField(record, 'OETraj', Types.Record)
            opt['Logger'].info('{} found OETraj record'.format(system_title))
            ligTraj = utl.RequestOEField( oetrajRecord, 'LigTraj', Types.Chem.Mol)
            opt['Logger'].info('{} #atoms, #confs in ligand traj OEMol: {}, {}'
                               .format(system_title, ligTraj.NumAtoms(), ligTraj.NumConfs()))
            protTraj = utl.RequestOEField( oetrajRecord, 'ProtTraj', Types.Chem.Mol)
            opt['Logger'].info('{} #atoms, #confs in protein traj OEMol: {}, {}'
                               .format(system_title, protTraj.NumAtoms(), protTraj.NumConfs()))

            # Compute PBSA energies for the protein-ligand complex
            zapBind, zapBindPB, zapDesolEl, zapIntEl, zapBindSA25, saBuried = mmpbsa.TrajPBSA(
                                       ligTraj, protTraj)
            if zapBind is None:
                raise ValueError('{} Calculation of PBSA energies failed'.format(system_title) )
            # generate Surface Areas energy for buried SA based on 0.006 kcal/mol/A^2
            zapBindSA6 = [sa*-0.006 for sa in saBuried]

            # Create new record with traj interaction energy results
            opt['Logger'].info('{} writing trajPBSA OERecord'.format(system_title) )
            trajPBSA = OERecord()
            #
            zapBind_field = OEField("OEZap_PBSA25_Bind", Types.FloatVec,
                                    meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajPBSA.set_value(zapBind_field, zapBind)
            #
            zapBindPB_field = OEField("OEZap_PB_Bind", Types.FloatVec,
                                      meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajPBSA.set_value(zapBindPB_field, zapBindPB)
            #
            zapDesolEl_field = OEField("OEZap_PB_Desolvation", Types.FloatVec,
                                       meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajPBSA.set_value(zapDesolEl_field, zapDesolEl)
            #
            zapIntEl_field = OEField("OEZap_PB_Interaction", Types.FloatVec,
                                     meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajPBSA.set_value(zapIntEl_field, zapIntEl)
            #
            zapBindSA25_field = OEField("OEZap_SA25_Bind", Types.FloatVec,
                                        meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajPBSA.set_value(zapBindSA25_field, zapBindSA25)
            #
            saBuried_field = OEField("OEZap_BuriedArea", Types.FloatVec)
            trajPBSA.set_value(saBuried_field, saBuried)
            #
            zapBindSA6_field = OEField("OEZap_SA6_Bind", Types.FloatVec,
                                       meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            trajPBSA.set_value(zapBindSA6_field, zapBindSA6)

            # If the OETraj Interaction Energies has been done calculate MMPBSA values
            if 'TrajIntE' in analysesDone:
                opt['Logger'].info('{} found TrajIntE analyses'.format(system_title) )

                # Extract the relevant P-L Interaction Energies from the record
                oeTrjIntERecord = utl.RequestOEField(record, 'TrajIntE', Types.Record)
                opt['Logger'].info('{} found TrajIntE record'.format(system_title))
                PLIntE = utl.RequestOEField(oeTrjIntERecord,
                                            'protein_ligand_interactionEnergy', Types.FloatVec)
                opt['Logger'].info('{} found Protein-Ligand force field interaction energies'
                                   .format(system_title))
                # Calculate  and store MMPB and MMPBSA energies on the trajPBSA record
                zapMMPB = [eInt+eDesol for eInt,eDesol in zip(PLIntE, zapDesolEl)]
                zapMMPB_field = OEField("OEZap_MMPB_Bind", Types.FloatVec,
                                        meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
                trajPBSA.set_value(zapMMPB_field, zapMMPB)
                #
                zapMMPBSA = [eMMPB+eSA6 for eMMPB,eSA6 in zip(zapMMPB, zapBindSA6)]
                zapMMPBSA_field = OEField("OEZap_MMPBSA6_Bind", Types.FloatVec,
                                          meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
                trajPBSA.set_value(zapMMPBSA_field, zapMMPBSA)

                # Average MMPBSA for all the trajectory frames
                np_mmpbsa = np.array(zapMMPBSA)
                avg_mmpbsa = np_mmpbsa.mean()
                std_mmpbsa = np_mmpbsa.std()

                # Add to the record the Average MMPBSA floe report label
                record.set_value(Fields.floe_report_label, "MMPBSA = {:.1f}  &plusmn; {:.1f} kcal/mol".
                                 format(avg_mmpbsa, std_mmpbsa))

            # Add the trajPBSA record to the parent record
            record.set_value(OEField('TrajPBSA', Types.Record), trajPBSA)
            analysesDone.append('TrajPBSA')
            record.set_value(OEField('AnalysesDone', Types.StringVec), analysesDone)
            opt['Logger'].info('{} finished writing TrajPBSA OERecord'.format(system_title))

            self.success.emit(record)

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} in TrajPBSACube'.format(str(e)))
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


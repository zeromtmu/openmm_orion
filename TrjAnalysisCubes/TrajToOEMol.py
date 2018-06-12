import os
import traceback

from floe.api import (ParallelMixin,
                      parameter)

# Just for old orion testing
from datarecord import OEField, Types, OERecord


from cuberecord import OERecordComputeCube
from Standards import (Fields,
                       MDRecords,
                       MDStageNames)

from openeye import oechem

import TrjAnalysisCubes.TrajToOEMol_utils as trjutl

def CheckAndGetValue( record, field, rType):
    if not record.has_value(OEField(field,rType)):
        #opt['Logger'].warn('Missing record field {}'.format( field))
        print( 'Missing record field {}'.format( field))
        raise ValueError('The record does not have field {}'.format( field))
    else:
        return record.get_value(OEField(field,rType))


class TrajToOEMolCube(ParallelMixin, OERecordComputeCube):
    title = 'Traj to OEMol Cube'

    version = "0.0.0"
    classification = [["Simulation", "Traj Analysis"]]
    tags = ['Parallel Cube']

    description = """
    Converting MD Traj into multiconf OEMols for Ligand and Protein

    This cube will take in the MD traj file containing
    the solvated protein:ligand complex and extract
    multiconf OEMols for Ligand and Protein.

    Input parameters:
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        return

    def process(self, record, port):
        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)

            # Logger string
            opt['Logger'].info(' ')
            system_title = CheckAndGetValue( record, 'Title_PLMD', Types.String)
            opt['Logger'].info('{}: Attempting MD Traj conversion into OEMols'
                .format(system_title) )
            floeID = CheckAndGetValue( record, 'ID_PLMD', Types.Int)
            opt['Logger'].info('{} floe ID: {}'.format(system_title, floeID) )

            # Extract the MDStageRecord list
            md_stages = CheckAndGetValue( record, 'MDStages', Types.RecordVec)
            if len(md_stages)<2:
                raise ValueError('{} does not have at least 2 MD Stages'.format(system_title) )

            # Extract and verify the traj filename for the last MD stage
            md_stageLast_record = md_stages[-1]
            lastName = CheckAndGetValue( md_stageLast_record, 'Stage_name', Types.String)
            if lastName!='NPT':
                raise ValueError('Cannot find the NPT stage')
            trajName = CheckAndGetValue( md_stageLast_record, 'Trajectory', Types.String)
            if trajName==None:
                opt['Logger'].info('Traj name was None' )
                raise ValueError('{} invalid Trajectory filename:'.format(system_title) )
            elif trajName=='':
                opt['Logger'].info('Traj name was an empty string' )
                raise ValueError('invalid Trajectory filename:' )
            opt['Logger'].info('{} Trajectory filename: {}'.format(system_title,trajName) )

            if os.path.exists( trajName):
                opt['Logger'].info('{} Trajectory file {} exists'.format(system_title,trajName) )
            else:
                opt['Logger'].info('{} cannot find Trajectory file {}'.format(system_title,trajName) )
                raise ValueError('{}: Cannot find Trajectory file'.format(system_title) )

            # Extract the Setup Topology
            md_stage0_record = md_stages[0]
            setupName = CheckAndGetValue( md_stage0_record, 'Stage_name', Types.String)
            if setupName!='SETUP':
                raise ValueError('Cannot find the SETUP stage')
            md_system = CheckAndGetValue( md_stage0_record, 'MDSystem', Types.Record)
            setupOEMol = CheckAndGetValue( md_system, 'Topology_OEMol', Types.Chem.Mol)
            opt['Logger'].info('Setup topology has {} atoms'.format(setupOEMol.NumAtoms()))

            # Generate multiconformer protein and ligand OEMols from the trajectory
            opt['Logger'].info('{} Generating protein and ligand trajectory OEMols'
                .format( system_title))
            ptraj, ltraj = trjutl.ExtractAlignedProtLigTraj_hdf5(setupOEMol, trajName )
            opt['Logger'].info('{} #atoms, #confs in protein traj OEMol: {}, {}'
                .format( system_title, ptraj.NumAtoms(), ptraj.NumConfs()) )
            opt['Logger'].info('{} #atoms, #confs in ligand traj OEMol: {}, {}'
                .format( system_title, ltraj.NumAtoms(), ltraj.NumConfs()) )

            # Generate average and median protein and ligand OEMols from ptraj, ltraj
            opt['Logger'].info('{} Generating protein and ligand median and average OEMols'
                .format( system_title))
            ligMedian, protMedian, ligAverage, protAverage = trjutl.AnalyseProteinLigandTrajectoryOEMols( ltraj, ptraj)

            # Overwrite MDStages with only first (setup) and last (production) stages
            newMDStages = [ md_stage0_record, md_stageLast_record]
            record.set_value( OEField( 'MDStages', Types.RecordVec), newMDStages)

            # Create new record with OETraj results
            oetrajRecord = OERecord()
            oetrajRecord.set_value( OEField( 'ProtTraj', Types.Chem.Mol), ptraj)
            oetrajRecord.set_value( OEField( 'LigTraj', Types.Chem.Mol), ltraj)
            oetrajRecord.set_value( OEField( 'LigMedian', Types.Chem.Mol), ligMedian)
            oetrajRecord.set_value( OEField( 'ProtMedian', Types.Chem.Mol), protMedian)
            oetrajRecord.set_value( OEField( 'LigAverage', Types.Chem.Mol), ligAverage)
            oetrajRecord.set_value( OEField( 'ProtAverage', Types.Chem.Mol), protAverage)
            record.set_value( OEField( 'OETraj', Types.Record), oetrajRecord)

            analysesDone = None
            try:
                analysesDone = CheckAndGetValue( record, 'AnalysesDone', Types.StringVec)
                opt['Logger'].info('{}: found AnalysesDone list'.format( system_title) )
                analysesDone.append( 'OETraj')
            except:
                analysesDone = [ 'OETraj' ]
                opt['Logger'].info('{}: created AnalysesDone list'.format( system_title) )
            if analysesDone is None:
                raise ValueError('{} AnalysesDone list does not exist'.format( system_title) )
            record.set_value( OEField( 'AnalysesDone', Types.StringVec), analysesDone)
            opt['Logger'].info('{}: saved protein and ligand traj OEMols'.format( system_title) )


            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


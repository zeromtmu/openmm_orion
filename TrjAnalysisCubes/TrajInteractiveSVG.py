import os
import traceback

from floe.api import (ParallelMixin,
                      parameter)

# Just for old orion testing
from datarecord import (Types,
                        Meta,
                        OEFieldMeta,
                        OEField,
                        OERecord)


from cuberecord import OERecordComputeCube
from Standards import (Fields,
                       MDRecords,
                       MDStageNames)

from openeye import oechem

import ensemble2img

def CheckAndGetValue( record, field, rType):
    if not record.has_value(OEField(field,rType)):
        #opt['Logger'].warn('Missing record field {}'.format( field))
        print( 'Missing record field {}'.format( field))
        raise ValueError('The record does not have field {}'.format( field))
    else:
        return record.get_value(OEField(field,rType))


class TrajInteractiveSVG(ParallelMixin, OERecordComputeCube):
    title = 'Generate Interactive SVG of Trajectory Analysis'

    version = "0.0.1"
    classification = [["Simulation", "Traj Analysis"]]
    tags = ['Parallel Cube']

    description = """
    Generate Interactive SVG of Trajectory Analysis

    This cube will take in the MD traj OEMols containing
    the protein and ligand components of the complex and generate
    an interactive svg showing protein-ligand interactions, calculated
    B-factors, and torsional fluctuations of the bound ligand.

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
            opt['Logger'].info('{} Attempting to generate interative svg of Traj'
                .format(system_title) )
            floeID = CheckAndGetValue( record, 'ID_PLMD', Types.Int)
            opt['Logger'].info('{} floe ID: {}'.format(system_title, floeID) )

            # Check that the OETraj analysis has been done
            analysesDone = CheckAndGetValue( record, 'AnalysesDone', Types.StringVec)
            if 'OETraj' not in analysesDone:
                raise ValueError('{} does not have OETraj analyses done'.format(system_title) )
            else:
                opt['Logger'].info('{} found OETraj analyses'.format(system_title) )

            # Extract the relevant traj OEMols from the OETraj record
            oetrajRecord = CheckAndGetValue( record, 'OETraj', Types.Record)
            opt['Logger'].info('{} found OETraj record'.format(system_title) )
            ligTraj = CheckAndGetValue( oetrajRecord, 'LigTraj', Types.Chem.Mol)
            opt['Logger'].info('{} #atoms, #confs in ligand traj OEMol: {}, {}'
                .format( system_title, ligTraj.NumAtoms(), ligTraj.NumConfs()) )
            protTraj = CheckAndGetValue( oetrajRecord, 'ProtTraj', Types.Chem.Mol)
            opt['Logger'].info('{} #atoms, #confs in protein traj OEMol: {}, {}'
                .format( system_title, protTraj.NumAtoms(), protTraj.NumConfs()) )
            ligAvg = CheckAndGetValue( oetrajRecord, 'LigAverage', Types.Chem.Mol)
            opt['Logger'].info('{} #atoms, #confs in ligand Average OEMol: {}, {}'
                .format( system_title, ligAvg.NumAtoms(), ligAvg.NumConfs()) )
            protAvg = CheckAndGetValue( oetrajRecord, 'ProtAverage', Types.Chem.Mol)
            opt['Logger'].info('{} #atoms, #confs in protein Average OEMol: {}, {}'
                .format( system_title, protAvg.NumAtoms(), protAvg.NumConfs()) )

            # Cluster ligand traj into cluster OEMols with matching protein OEMols
            opt['Logger'].info('{} Generating svg'.format(system_title) )
            trajSVGsvg = ensemble2img.run_ensemble2img( ligAvg, protAvg, ligTraj, protTraj)

            # Create new record with trajSVG results
            opt['Logger'].info('{} writing trajSVG OERecord'.format(system_title) )
            trajSVG = OERecord()
            trajSVG_field = OEField( 'TrajSVG', Types.String, meta=OEFieldMeta().set_option(Meta.Hints.Image_SVG))
            trajSVG.set_value( trajSVG_field, trajSVGsvg)

            #
            record.set_value( OEField( 'TrajSVG', Types.Record), trajSVG)
            analysesDone.append( 'TrajSVG')
            record.set_value( OEField( 'AnalysesDone', Types.StringVec), analysesDone)
            opt['Logger'].info('{} finished writing trajSVG OERecord'.format(system_title) )

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


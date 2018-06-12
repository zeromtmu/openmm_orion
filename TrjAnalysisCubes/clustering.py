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

import TrjAnalysisCubes.clustering_utils as clusutl

def CheckAndGetValue( record, field, rType):
    if not record.has_value(OEField(field,rType)):
        #opt['Logger'].warn('Missing record field {}'.format( field))
        print( 'Missing record field {}'.format( field))
        raise ValueError('The record does not have field {}'.format( field))
    else:
        return record.get_value(OEField(field,rType))


class ClusterOETrajCube(ParallelMixin, OERecordComputeCube):
    title = 'Cluster Protein-Ligand Traj OEMols'

    version = "0.0.1"
    classification = [["Simulation", "Traj Analysis"]]
    tags = ['Parallel Cube']

    description = """
    Cluster  multiconf MD trajectory OEMols for Ligand and Protein

    This cube will take in the MD traj OEMols containing
    the protein and ligand components of the complex and cluster
    them based on ligand RMSD.

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
            opt['Logger'].info('{}: Attempting to cluster MD Traj conversion protein-ligand OEMols'
                .format(system_title) )
            floeID = CheckAndGetValue( record, 'ID_PLMD', Types.Int)
            opt['Logger'].info('{} floe ID: {}'.format(system_title, floeID) )

            # Extract the MDStageRecord list
            analysesDone = CheckAndGetValue( record, 'AnalysesDone', Types.StringVec)
            if 'OETraj' not in analysesDone:
                raise ValueError('{} does not have OETraj analyses done'.format(system_title) )
            else:
                opt['Logger'].info('{} found OETraj analyses'.format(system_title) )

            # Extract and verify the traj filename for the last MD stage

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


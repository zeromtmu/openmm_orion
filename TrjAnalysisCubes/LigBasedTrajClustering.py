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

#import TrjAnalysisCubes.Clustering_utils as clusutl
import oetrajanalysis.Clustering_utils as clusutl

def CheckAndGetValueFull( record, field, rType):
    if not record.has_value(OEField(field,rType)):
        #opt['Logger'].warn('Missing record field {}'.format( field))
        print( 'Missing record field {}'.format( field))
        raise ValueError('The record does not have field {}'.format( field))
    else:
        return record.get_value(OEField(field,rType))

def CheckAndGetValue( record, field):
    if not record.has_value(field):
        #opt['Logger'].warn('Missing record field {}'.format( field))
        print( 'Missing record field {}'.format( field.get_name() ))
        raise ValueError('The record does not have field {}'.format( field.get_name() ))
    else:
        return record.get_value(field)


class ClusterOETrajCube(ParallelMixin, OERecordComputeCube):
    title = 'Cluster Protein-Ligand Traj OEMols'

    version = "0.1.0"
    classification = [["Simulation", "Traj Analysis"]]
    tags = ['Parallel Cube']

    description = """
    Cluster  multiconf MD trajectory OEMols for Ligand and Protein

    This cube will take in the MD traj OEMols containing
    the protein and ligand components of the complex and cluster
    them based on ligand RMSD.

    Input parameters:
    -------
    oechem.OEDataRecord - Streamed-in input data for the system to cluster

    Output:
    -------
    oechem.OEDataRecord - Stream of output data for the clustered system
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
            opt['Logger'].info(' Beginning ClusterOETrajCube')
            system_title = CheckAndGetValue( record, Fields.title)
            opt['Logger'].info('{} Attempting to cluster MD Traj'
                .format(system_title) )

            # Check that the OETraj analysis has been done
            analysesDone = CheckAndGetValueFull( record, 'AnalysesDone', Types.StringVec)
            if 'OETraj' not in analysesDone:
                raise ValueError('{} does not have OETraj analyses done'.format(system_title) )
            else:
                opt['Logger'].info('{} found OETraj analyses'.format(system_title) )

            # Extract the relevant traj OEMols from the OETraj record
            oetrajRecord = CheckAndGetValueFull( record, 'OETraj', Types.Record)
            opt['Logger'].info('{} found OETraj record'.format(system_title) )
            ligTraj = CheckAndGetValueFull( oetrajRecord, 'LigTraj', Types.Chem.Mol)
            opt['Logger'].info('{} #atoms, #confs in ligand traj OEMol: {}, {}'
                .format( system_title, ligTraj.NumAtoms(), ligTraj.NumConfs()) )

            # Cluster ligand traj into cluster OEMols with matching protein OEMols
            opt['Logger'].info('{} starting clustering'.format(system_title) )
            clusResults = clusutl.ClusterLigTraj( ligTraj)
            opt['Logger'].info('{} plotting RMSD histogram'.format(system_title) )
            trajHistRMSD_svg = clusutl.ClusterLigTrajHistRMSD(clusResults)
            # Heat map is large, time-consuming, and non-essential... put back if needed
            #opt['Logger'].info('{} plotting RMSD heat map'.format(system_title) )
            #trajHeatRMSD_png = clusutl.ClusterLigTrajHeatRMSD(clusResults)
            opt['Logger'].info('{} plotting cluster strip plot'.format(system_title) )
            trajClus_svg = clusutl.ClusterLigTrajClusPlot(clusResults)

            # Create new record with trajClus results
            opt['Logger'].info('{} writing trajClus OERecord'.format(system_title) )
            trajClus = OERecord()
            #
            ClusterMethod_field = OEField( 'ClusterMethod', Types.String)
            trajClus.set_value( ClusterMethod_field, clusResults['ClusterMethod'])
            #
            HDBSCAN_alpha_field = OEField( 'HDBSCAN_alpha', Types.Float)
            trajClus.set_value( HDBSCAN_alpha_field, clusResults['HDBSCAN_alpha'])
            #
            nFrames_field = OEField( 'nFrames', Types.Int)
            trajClus.set_value( nFrames_field, clusResults['nFrames'])
            #
            nClusters_field = OEField( 'nClusters', Types.Int)
            trajClus.set_value( nClusters_field, clusResults['nClusters'])
            #
            clusterCounts_field = OEField( 'ClusterCounts', Types.IntVec)
            trajClus.set_value( clusterCounts_field, clusResults['ClusterCounts'])
            #
            Clusters_field = OEField( 'Clusters', Types.IntVec)
            trajClus.set_value( Clusters_field, clusResults['ClusterVec'])
            #
            HistSVG_field = OEField( 'HistSVG', Types.String, meta=OEFieldMeta().set_option(Meta.Hints.Image_SVG))
            trajClus.set_value( HistSVG_field, trajHistRMSD_svg)
            #
            ClusSVG_field = OEField( 'ClusSVG', Types.String, meta=OEFieldMeta().set_option(Meta.Hints.Image_SVG))
            trajClus.set_value( ClusSVG_field, trajClus_svg)
            # Heat map is large, time-consuming, and non-essential... put back if needed
            # HeatPNG_field = OEField( 'HeatPNG', Types.Blob, meta=OEFieldMeta().set_option(Meta.Hints.Image_PNG))
            # trajClus.set_value( HeatPNG_field, trajHeatRMSD_png)

            #
            record.set_value( OEField( 'TrajClus', Types.Record), trajClus)
            analysesDone.append( 'TrajClus')
            record.set_value( OEField( 'AnalysesDone', Types.StringVec), analysesDone)
            opt['Logger'].info('{} finished writing trajClus OERecord'.format(system_title) )

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


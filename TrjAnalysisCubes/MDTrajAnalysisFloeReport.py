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


def CheckAndGetValue( record, field, rType):
    if not record.has_value(OEField(field,rType)):
        #opt['Logger'].warn('Missing record field {}'.format( field))
        print( 'Missing record field {}'.format( field))
        raise ValueError('The record does not have field {}'.format( field))
    else:
        return record.get_value(OEField(field,rType))


class MDTrajAnalysisToFiles(ParallelMixin, OERecordComputeCube):
    title = 'Extract relevant outputs of MD Traj Analysis'

    version = "0.0.1"
    classification = [["Simulation", "Traj Analysis"]]
    tags = ['Parallel Cube']

    description = """
    Extract relevant outputs of Ligand and Protein
    Short Traj MD Traj Analysis and write them to files.

    This cube takes as input the OERecord containing the work
    product of trajectory analysis on Short Traj MD results.
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

            # title of entire solvated protein-ligand system
            opt['Logger'].info(' ')
            system_title = CheckAndGetValue( record, 'Title_PLMD', Types.String)
            opt['Logger'].info('{} Attempting to extract MD Traj Analysis results'
                .format(system_title) )
            floeID = CheckAndGetValue( record, 'ID_PLMD', Types.Int)
            opt['Logger'].info('{} floe ID: {}'.format(system_title, floeID) )
            ligInitPose = CheckAndGetValue( record, 'Ligand', Types.Chem.Mol)

            # Extract the traj SVG from the OETraj record
            analysesDone = CheckAndGetValue( record, 'AnalysesDone', Types.StringVec)
            if 'OETraj' not in analysesDone:
                raise ValueError('{} does not have OETraj analyses done'.format(system_title) )
            else:
                opt['Logger'].info('{} found OETraj analyses'.format(system_title) )
            # Extract the relevant traj SVG from the OETraj record
            oetrajRecord = CheckAndGetValue( record, 'OETraj', Types.Record)
            opt['Logger'].info('{} found OETraj record'.format(system_title) )
            trajSVG = CheckAndGetValue( oetrajRecord, 'TrajSVG', Types.String)

            # Extract the three plots from the TrajClus record
            analysesDone = CheckAndGetValue( record, 'AnalysesDone', Types.StringVec)
            if 'TrajClus' not in analysesDone:
                raise ValueError('{} does not have TrajClus analyses done'.format(system_title) )
            else:
                opt['Logger'].info('{} found TrajClus analyses'.format(system_title) )
            # Extract the relevant traj SVG from the TrajClus record
            clusRecord = CheckAndGetValue( record, 'TrajClus', Types.Record)
            opt['Logger'].info('{} found TrajClus record'.format(system_title) )
            trajHistRMSD_svg = CheckAndGetValue( clusRecord, 'HistSVG', Types.String)
            trajClus_svg = CheckAndGetValue( clusRecord, 'ClusSVG', Types.String)
            #trajHeatRMSD_png = CheckAndGetValue( clusRecord, 'HeatPNG', Types.Blob)
            opt['Logger'].info('{} found the TrajClus plots'.format(system_title) )

            # write files for each result
            with oechem.oemolostream(system_title+'_ligInitPose.oeb') as ofs:
                oechem.OEWriteConstMolecule(ofs,ligInitPose)
            with open(system_title+'_traj.svg','w') as ofs:
                ofs.write( trajSVG)
            with open(system_title+'_histRMSD.svg','w') as ofs:
                ofs.write( trajHistRMSD_svg)
            #with open(system_title+'_heatRMSD.png','wb') as ofs:
            #    ofs.write( trajHeatRMSD_png)
            with open(system_title+'_clusters.svg','w') as ofs:
                ofs.write( trajClus_svg)

            # write useful text to go into floe report
            with open(system_title+'_trajAnalysis.txt','w') as ofs:
                ofs.write( '{} : Analysis of Short Trajectory MD\n'.format(ligInitPose.GetTitle()) )
                ofs.write( '{} : has {} atoms\n'.
                    format( ligInitPose.GetTitle(), ligInitPose.NumAtoms() ) )
                nFrames = CheckAndGetValue( clusRecord, 'nFrames', Types.Int)
                ofs.write( 'Clustering ligand trajectory of {} frames\n'.format( nFrames))
                ofs.write( '    based on active site alignment:\n')
                clusMethod = CheckAndGetValue( clusRecord, 'ClusterMethod', Types.String)
                alpha = CheckAndGetValue( clusRecord, 'HDBSCAN_alpha', Types.Float)
                ofs.write( 'Using clustering method {} with alpha {}\n'.format( clusMethod, alpha))
                nClusters = CheckAndGetValue( clusRecord, 'nClusters', Types.Int)
                ofs.write( 'produced {} clusters\n'.format( nClusters))
                clusCounts = CheckAndGetValue( clusRecord, 'ClusterCounts', Types.IntVec)
                for i, count in enumerate(clusCounts):
                    ofs.write( 'cluster {} contains {} frames'.format( i, count))

            opt['Logger'].info('{} finished writing analysis files'.format(system_title) )

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


# (C) 2019 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

from floe.api.orion import in_orion

from MDOrion.Standards.utils import ParmedData, MDStateData

from datarecord import OEPrimaryMolField

from datarecord import (Types,
                        Meta,
                        OEFieldMeta,
                        OEField,
                        OERecord)

import os

from MDOrion.Standards import utils

# ------------ Stage Standard Names ------------- #

class MDStageTypes:
    SETUP = 'SETUP'
    MINIMIZATION = 'MINIMIZATION'
    NVT = 'NVT'
    NPT = 'NPT'
    FEC = 'FEC'


# ------------ MD Engines ------------- #

class MDEngines:
    OpenMM = 'OpenMM'
    Gromacs = 'Gromacs'
    all = [OpenMM, Gromacs]


# ---------------- File  Name Standards -------------- #

class MDFileNames:
    topology = 'topology.oeb'
    state = 'state.pickle'
    trajectory = "trajectory.tar.gz"
    trajectory_conformers = "trajectory.oeb"
    mddata = "data.tar.gz"

# ---------------- Field Standards -------------- #


class Fields:

    # The Title field is used to set the system name
    title = OEField("Title_OPLMD", Types.String, meta=OEFieldMeta().set_option(Meta.Source.ID))

    # The ID field should be used as identification number for ligands, proteins or complexes
    id = OEField("ID_OPLMD", Types.Int, meta=OEFieldMeta().set_option(Meta.Source.ID))

    # The SysID field is used to keep track of the system order
    sysid = OEField("SysID_OPLMD", Types.Int, meta=OEFieldMeta().set_option(Meta.Source.ID))

    # The ConfID field is used to identify a particular conformer
    confid = OEField("ConfID_OPLMD", Types.Int, meta=OEFieldMeta().set_option(Meta.Source.ID))

    # The Ligand field should be used to save in a record a ligand as an OEMolecule
    ligand = OEField("Ligand_OPLMD", Types.Chem.Mol, meta=OEFieldMeta().set_option(Meta.Hints.Chem.Ligand))

    # The ligand name
    ligand_name = OEField("Ligand_name_OPLMD", Types.String)

    # The Protein field should be used to save in a record a Protein as an OEMolecule
    protein = OEField("Protein_OPLMD", Types.Chem.Mol, meta=OEFieldMeta().set_option(Meta.Hints.Chem.Protein))

    # Primary Molecule
    primary_molecule = OEPrimaryMolField()

    # Parmed Structure Field
    pmd_structure = OEField('Structure_Parmed_OPLMD', ParmedData)

    # The Stage Name
    stage_name = OEField('Stage_name_OPLMD', Types.String)

    # The Stage Type
    stage_type = OEField('Stage_type_OPLMD', Types.String)

    # Topology Field
    topology = OEField('Topology_OPLMD', Types.Chem.Mol, meta=OEFieldMeta().set_option(Meta.Hints.Chem.PrimaryMol))

    # Log Info
    log_data = OEField('Log_data_OPLMD', Types.String)

    # MD State
    md_state = OEField("MDState_OPLMD", MDStateData)

    # # MD System Field
    # md_system = OEField("MDSystem_OPLMD", Types.Record)

    # Trajectory
    if in_orion():
        trajectory = OEField("Trajectory_OPLMD", Types.Int)
    else:
        trajectory = OEField("Trajectory_OPLMD", Types.String)

    # MD Data
    if in_orion():
        mddata = OEField("MDData_OPLMD", Types.Int)
    else:
        mddata = OEField("MDData_OPLMD", Types.String)

    # This Field is introduced to deal with record trajectory field
    # that are linked to Orion S3 storage but they are locally used
    # orion_local_trj_field = OEField("Trajectory_OPLMD", Types.Int)

    # Collection is used to offload data from the record which mush be < 100Mb
    collection = OEField("Collection_OPLMD", Types.String)

    # Stage list Field
    md_stages = OEField("MDStages_OPLMD", Types.RecordVec)

    # Analysis Fields
    free_energy = OEField('FE_OPLMD', Types.Float,
                          meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))

    free_energy_err = OEField('FE_Error_OPLMD', Types.Float,
                              meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))

    floe_report = OEField('Floe_report_OPLMD', Types.String)

    floe_report_svg_lig_depiction = OEField("Floe_report_lig_svg_OPLMD", Types.String,
                                            meta=OEFieldMeta().set_option(Meta.Hints.Image_SVG))

    floe_report_label = OEField('Floe_report_label_OPLMD', Types.String)



# ---------------- Record Standards -------------- #

class MDRecords:
    @staticmethod
    def MDSystemRecord(molecule, state):
        record = OERecord()
        record.set_value(Fields.topology, molecule)
        record.set_value(Fields.md_state, state)
        return record

    @staticmethod
    def MDStageRecord(stage_name,
                      stage_type,
                      system_record,
                      log=None,
                      trajectory=None,
                      trajectory_engine=None,
                      orion_name="OrionFile"):

        record = OERecord()

        record.set_value(Fields.stage_name, stage_name)
        record.set_value(Fields.stage_type, stage_type)
        record.set_value(Fields.md_system, system_record)

        if log is not None:
            record.set_value(Fields.log_data, log)
        if trajectory is not None:

            if trajectory_engine not in MDEngines.all:
                raise ValueError("The selected MD engine is not supported")

            trj_meta = OEFieldMeta()
            trj_meta.set_attribute(Meta.Annotation.Description, trajectory_engine)

            trj_field = OEField(Fields.trajectory.get_name(),
                                Fields.trajectory.get_type(),
                                meta=trj_meta)

            if not os.path.isfile(trajectory):
                raise IOError("The trajectory file has not been found: {}".format(trajectory))

            lf = utils.upload_file(trajectory, orion_name=orion_name)

            record.set_value(trj_field, lf)

        return record

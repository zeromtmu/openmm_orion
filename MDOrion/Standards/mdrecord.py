from MDOrion.Standards import Fields, MDRecords, MDEngines

from MDOrion.Standards import utils

import parmed

import copy

import os

import tarfile

from datarecord import (Meta,
                        OEFieldMeta,
                        OEField)
import datarecord

from openeye import oechem


def mdstages(f):
    def wrapper(*pos, **named):
        rec = pos[0]
        if not rec.has_field(Fields.md_stages):
            raise ValueError("The MD record does not have MD stages")

        return f(*pos, **named)

    return wrapper


class MDDataRecord(object):
    def __init__(self, record, inplace=True):
        if inplace:
            self.rec = record
        else:
            self.rec = copy.deepcopy(record)

    @property
    def get_record(self):
        """
        This method returns the record

        Return:
        -------
        record: OERecord
            The record
        """
        return self.rec

    @property
    def get_primary(self):
        """
        This method returns the primary molecule present on the record

        Return:
        -------
        record: OEMol
            The Primary Molecule
        """

        if not self.rec.has_field(Fields.primary_molecule):
            raise ValueError("The Primary Molecule has not been found on the record")

        return self.rec.get_value(Fields.primary_molecule)

    def set_primary(self, primary_mol):
        """
        This method sets the primary molecule on the record

        Return:
        -------
        record: Bool
            True if the primary molecule has been set on the record
        """

        if not isinstance(primary_mol, oechem.OEMol):
            raise ValueError("The Primary Molecule is not a valid OEMol: {}".format(primary_mol))

        self.rec.set_value(Fields.primary_molecule, primary_mol)

        return True

    @property
    def get_id(self):
        """
        This method returns ID field present on the record

        Return:
        -------
            : Int
            The record ID
        """

        if not self.rec.has_field(Fields.id):
            raise ValueError("The ID Field has not been found on the record")

        return self.rec.get_value(Fields.id)

    @property
    def has_id(self):
        """
        This method checks if the ID field is present on the record

        Return:
        -------
            : Bool
            True if the ID field is resent on the record otherwise False
        """

        if not self.rec.has_field(Fields.id):
            return False
        else:
            return True

    def set_id(self, id):
        """
        This method sets the ID field on the record

        Return:
        -------
            : Bool
            True if the ID has been set on the record
        """

        if not isinstance(id, int):
            raise ValueError(" The id must be an integer: {}".format(id))

        self.rec.set_value(Fields.id, id)

        return True

    @property
    def get_title(self):
        """
        This method returns the system title present on the record

        Return:
        -------
            : String
            The system title tring
        """

        if not self.rec.has_field(Fields.title):
            raise ValueError("The Title Field has not been found on the record")

        return self.rec.get_value(Fields.title)

    @property
    def has_title(self):
        """
        This method checks if the Title field is present on the record

        Return:
        -------
            : Bool
            True if the Title field is resent on the record otherwise False
        """

        if not self.rec.has_field(Fields.title):
            return False
        else:
            return True

    def set_title(self, title):
        """
        This method sets the system Title field on the record

        Return:
        -------
            : Bool
            True if the system Title has been set on the record
        """

        if not isinstance(title, str):
            raise ValueError(" The title must be a sting: {}".format(title))

        self.rec.set_value(Fields.title, title)

        return True

    @property
    @mdstages
    def get_last_stage(self):
        """
        This method returns the last MD stage of the MD record stages

        Return:
        -------
        record: OERecord
            The last stage of the MD record stages
        """

        return self.rec.get_value(Fields.md_stages)[-1]

    @mdstages
    def get_stage_by_name(self, name='last'):
        """
        This method returns a MD stage selected by passing the string stage name. If the
        string "last" is passed (default) the last MD stage is returned. If no stage name
        has been found an exception is raised.

        Parameters:
        -----------
        name: String
            The MD stage name

        Return:
        -------
        record: OERecord
            The MD stage selected by its name
        """
        stages = self.rec.get_value(Fields.md_stages)

        stg_names = []

        if name == 'last':
            return stages[-1]

        for stage in stages:
            stg_name = stage.get_value(Fields.stage_name)
            if stg_name == name:
                return stage
            else:
                stg_names.append(stg_name)

        raise ValueError("The Stage name has not been found: {} available names: {}".format(name, stg_names))

    @mdstages
    def delete_stage_by_name(self, name='last'):
        """
        This method deletes an MD stage selected by passing the string name. If the
        string "last" is passed (default) the last MD stage is deleted. If no stage name
        has been found an exception is raised.

        Parameters:
        -----------
        name: String
            The MD stage name

        Return:
        -------
        record: Bool
            True if the deletion was successful
        """
        stages = self.rec.get_value(Fields.md_stages)

        stg_names = []

        if name == 'last':
            del stages[-1]

        else:
            for idx in range(0, len(stages)):

                stg_name = stages[idx].get_value(Fields.stage_name)

                if stg_name == name:
                    del stages[idx]
                    self.rec.set_value(Fields.md_stages, stages)
                    return True
                else:
                    stg_names.append(stg_name)

        raise ValueError("The Stage name has not been found: {} available names: {}".format(name, stg_names))

    @mdstages
    def has_stage_name(self, name):
        """
        This method returns True if MD stage selected by passing the string name is present
        on the MD stage record otherwise False.

        Parameters:
        -----------
        name: String
            The MD stage name

        Return:
        -------
        record: Bool
            True if the MD stage name is present on the MD stages record otherwise False
        """

        stages = self.rec.get_value(Fields.md_stages)

        for stage in stages:
            stg_name = stage.get_value(Fields.stage_name)
            if stg_name == name:
                return True

        return False

    @mdstages
    def set_last_stage(self, mdstage):
        """
        This method overwrites the last MD stage presents in the MD stage record with the
        new passed md stage. If the setting was successful True is returned. If the passed MD
        stage is not valid an exception is raised.

        Parameters:
        -----------
        mdstage: OERecord
            The MD stage to be written on the MD stage record

        Return:
        -------
            : Bool
            True if the setting was successful
        """
        if not isinstance(mdstage, datarecord.datarecord.OERecord):
            raise ValueError("The passed md stage object is not a valid MD Stage: {}".format(type(mdstage)))

        stages = self.rec.get_value(Fields.md_stages)

        stages[-1] = mdstage

        self.rec.set_value(Fields.md_stages, stages)

        return True

    @mdstages
    def append_stage(self, mdstage):
        """
        This method appends a new MD stage on the MD stage record. If the appending
        was successful True is returned. If the passed MD stage is not valid an exception
        is raised.

        Parameters:
        -----------
        mdstage: OERecord
            The MD stage to be written on the MD stage record

        Return:
        -------
            : Bool
            True if the appending was successful
        """

        if not isinstance(mdstage, datarecord.datarecord.OERecord):
            raise ValueError("The passed md stage object is not a valid MD Stage: {}".format(type(mdstage)))

        stages = self.rec.get_value(Fields.md_stages)

        stages.append(mdstage)

        self.rec.set_value(Fields.md_stages, stages)

        return True

    @mdstages
    def get_stage_state(self, name='last', stage=None):
        """
        This method returns the MD State of the selected MD stage. If no stage name is passed
        the last stage is selected. If a MD stage is passed the MD state of the passed MD stage is returned.
        If a name and a MD stage is specified the MD stage name is ignored

        Parameters:
        -----------
        name: String
            The MD stage name
        stage: OERecord
            The MD stage record

        Return:
        -------
            : OERecord
            The MD state of the selected MD stage
        """

        if stage is not None:
            if not isinstance(stage, datarecord.datarecord.OERecord):
                raise ValueError("The passed md stage object is not a valid MD Stage: {}".format(type(stage)))
        else:
            stage = self.get_stage_by_name(name)

        md_system = stage.get_value(Fields.md_system)

        return md_system.get_value(Fields.md_state)

    @mdstages
    def get_stage_topology(self, name='last', stage=None):
        """
        This method returns the MD topology of the selected MD stage. If no stage name is passed
        the last stage is selected. If a MD stage is passed the MD topology of the passed MD stage is returned.
        If a name and a MD stage is specified the MD stage name is ignored

        Parameters:
        -----------
        name: String
            The MD stage name
        stage: OERecord
            The MD stage record

        Return:
        -------
            : OEMol
            The topology of the selected MD stage
        """

        if stage is not None:
            if not isinstance(stage, datarecord.datarecord.OERecord):
                raise ValueError("The passed md stage object is not a valid MD Stage: {}".format(type(stage)))
        else:
            stage = self.get_stage_by_name(name)

        md_system = stage.get_value(Fields.md_system)

        return md_system.get_value(Fields.topology)

    @mdstages
    def get_stage_info(self, name='last', stage=None):
        """
        This method returns the info of the selected MD stage. If no stage name is passed
        the last stage is selected. If a MD stage is passed the info of the passed MD stage are returned.
        If a name and a MD stage is specified the MD stage name is ignored

        Parameters:
        -----------
        name: String
            The MD stage name
        stage: OERecord
            The MD stage record

        Return:
        -------
            : String
            The info associated with the selected MD stage
        """

        if stage is not None:
            if not isinstance(stage, datarecord.datarecord.OERecord):
                raise ValueError("The passed md stage object is not a valid MD Stage: {}".format(type(stage)))
        else:
            stage = self.get_stage_by_name(name)

        return stage.get_value(Fields.log_data)

    @staticmethod
    def create_stage(stage_name,
                     stage_type,
                     topology,
                     mdstate,
                     log=None,
                     trajectory=None,
                     trajectory_engine=None,
                     orion_name="OrionFile"):
        """
        This method creates a new MD stage to be added to a MD stage record. If a trajectory file name
        is passed, the trajectory will be uploaded in AWS S3 in Orion

        Parameters:
        -----------
        stage_name: String
            The new MD stage name
        stage_type: String
            The MD stage type e.g. SETUP, MINIMIZATION etc.
        topology: OEMol
            The topology
        mdstate: OERecord
            The new mdstate made of state positions, velocities and box vectors
        log: String or None
            Log info
        trajectory: String, Int or None
            The trajectory name or id associated with the new MD stage
        trajectory_engine: String or None
            The MD engine used to generate the new MD stage. Possible names: OpenMM or Gromacs
        orion_name: String
            The Orion name to be visualize in the Orion UI

        Return:
        -------
            : OERecord
            The new created MD stage
        """

        md_system = MDRecords.MDSystemRecord(topology, mdstate)

        md_stage = MDRecords.MDStageRecord(stage_name,
                                           stage_type,
                                           md_system,
                                           log=log,
                                           trajectory=trajectory,
                                           trajectory_engine=trajectory_engine,
                                           orion_name=orion_name)
        return md_stage

    @mdstages
    def get_stage_trajectory(self, trajectory_name="trajectory.tar.gz", out_directory='.', name='last', stage=None):
        """
        This method returns the trajectory ID if the recovery of a trajectory linked to a MD stage is successful.
        If the MD record was locally generated the method returns the local trajectory file name if exists,
        otherwise an exception is raised. In this case trajectory_name and out_directory are ignored. If the
        MD record is generated in Orion and used in Orion or locally the method returns the trajectory ID if
        the trajectory can be downloaded from AWS S3 and the trajectory file is named as the passed parameter
        trajectory_name and downloaded in the passed out_directory. The stage to get the trajectory can be selected
        by its name by using the name parameter. If no name is specified the "last" stage is selected. If a MD stage
        is passed the trajectory associated with this stage is processed. If a name and a MD stage is specified
        the MD stage name is ignored

        Parameters:
        -----------
        trajectory_name: String
            The trajectory file name to be used in Orion or locally if the MD record was generated in Orion
        out_directory: String
            The output directory where saving the trajectory file in Orion or locally if the MD record was
            generated in Orion
        name: String
            The MD stage name
        stage: OERecord
            The MD stage record

        Return:
        -------
            : String or int
            Trajectory ID if the process was successful
        """

        if stage is not None:
            if not isinstance(stage, datarecord.datarecord.OERecord):
                raise ValueError("The passed md stage object is not a valid MD Stage: {}".format(type(stage)))
        else:
            stage = self.get_stage_by_name(name)

        if stage.get_value(Fields.trajectory) is None:
            return False

        if stage.has_value(Fields.trajectory):
            trj_id = stage.get_value(Fields.trajectory)

        elif stage.has_value(Fields.orion_local_trj_field):
            trj_id = stage.get_value(Fields.orion_local_trj_field)

        else:
            raise ValueError("No trajectory file have been found in the selected MD Stage {}".format(stage))

        trj_fn = os.path.join(out_directory, trajectory_name)

        fn_local = utils.download_file(trj_id, trj_fn, orion_delete=False)

        # with tarfile.open(fn_local) as tar:
        #     tar.extractall(path=out_directory)

        return fn_local

    @mdstages
    def set_stage_trajectory(self, trajectory_file_name, orion_name, trajectory_engine, name='last', stage=None):

        """
        This method returns True if a trajectory can be linked to a MD stage. If the trajectory file name
        linked to a trajectory file does not exist an exception is raised. If the MD record was locally generated
        the orion_name is ignored. If the trajectory is set in Orion the trajectory file is uploaded to S3 and an
        unique id related to the stack is set on the MD stage. The trajectory engine can be set as OpenMM or Gromacs.
        The stage to link the trajectory can be selected by its name by using the name parameter. If no name
        is specified the "last" stage is selected. If a MD stage is passed the trajectory associated with
        this stage is processed. If a name and a MD stage is specified the MD stage name is ignored.

        Parameters:
        -----------
        trajectory_file_name: String
            The trajectory file name to be used
        orion_name: String
            The Orion name to be visualize in the Orion UI
        name: String
            The MD stage name to link the trajectory
        stage: OERecord
            The MD stage record to link the trajectory

        Return:
        -------
            : Bool
            True if the trajectory was successful processed
        """

        if not os.path.isfile(trajectory_file_name):
            raise IOError("The trajectory file has not been found: {}".format(trajectory_file_name))

        lf = utils.upload_file(trajectory_file_name, orion_name=orion_name)

        if trajectory_engine not in MDEngines.all:
            raise ValueError("The selected MD engine is not supported")

        if trajectory_engine not in MDEngines.all:
            raise ValueError("The selected MD Engine {} is not supported. Valid are: {}".format(trajectory_engine,
                                                                                                MDEngines.all))
        trj_meta = OEFieldMeta()
        trj_meta.set_attribute(Meta.Annotation.Description, trajectory_engine)
        trj_field = OEField(Fields.trajectory.get_name(), Fields.trajectory.get_type(),
                            meta=trj_meta)

        if stage is not None:
            if not isinstance(stage, datarecord.datarecord.OERecord):
                raise ValueError("The passed md stage object is not a valid MD Stage: {}".format(type(stage)))

            stage.set_value(trj_field, lf)
            return True

        else:
            stage = self.get_stage_by_name(name)
            stage.set_value(trj_field, lf)

        stages = self.get_stages

        if name == 'last':
            stages[-1] = stage
        else:
            for idx in range(0, len(stages)):
                if stages[idx].get_value(Fields.stage_name) == name:
                    stages[idx] = stage
                    break

        self.rec.set_value(Fields.md_stages, stages)

        return True

    def init_stages(self, stage):
        """
        This method is used to initialize the record with a stage.

        Return:
        -------
            : Bool
            If the initialization succeeded True is returned
        """

        if not isinstance(stage, datarecord.datarecord.OERecord):
            raise ValueError("The passed md stage object is not a valid MD Stage: {}".format(type(stage)))

        self.rec.set_value(Fields.md_stages, [stage])

        return True

    @property
    @mdstages
    def get_stages(self):
        """
        This method returns the MD stage record list with all the MD stages.

        Return:
        -------
            : list
            The list of the MD stages
        """
        return self.rec.get_value(Fields.md_stages)

    @property
    @mdstages
    def get_stages_names(self):
        """
        This method returns the name list of the MD stages.

        Return:
        -------
            : list
            The list of the MD stage names
        """

        stages = self.rec.get_value(Fields.md_stages)
        stg_names = [stage.get_value(Fields.stage_name) for stage in stages]

        return stg_names

    @property
    def has_stages(self):
        """
        This method returns True if the record has a MD record list, False otherwise.

        Return:
        -------
            : Bool
            True if the record has a list of MD stages, False otherwise
        """
        if not self.rec.has_field(Fields.md_stages):
            return False
        else:
            return True

    @property
    def get_parmed(self):
        """
        This method returns the Parmed object present on the record. An exception is raised if
        the Parmed object cannot be found.

        Return:
        -------
            : Parmed Structure object
            The Parmed Structure object
        """

        if not self.rec.has_field(Fields.pmd_structure):
            raise ValueError("The Parmed Structure is not present on the record")

        return self.rec.get_value(Fields.pmd_structure)

    def set_parmed(self, pmdobj):
        """
        This method sets the Parmed object on the record. Return True if the setting was successful.

        Parameters:
        -----------
        pmjobj: Parmed Structure object
            The Parmed Structure object to be set on the record


        Return:
        -------
            : Bool
            True if the setting was successful
        """

        if not isinstance(pmdobj, parmed.Structure):
            raise ValueError("The passed Parmed object is not a valid Parmed Structure: {}".format(type(pmdobj)))

        self.rec.set_value(Fields.pmd_structure, pmdobj)

        return True

    @property
    def has_parmed(self):
        """
        This method checks if the Parmed object is on the record.

        Return:
        -------
            : Bool
            True if the Parmed object is on the record otherwise False
        """

        if not self.rec.has_field(Fields.pmd_structure):
            return False
        else:
            return True

    @property
    def delete_parmed(self):
        """
        This method deletes the Parmed object from the record. True is returned if the deletion was
        successful.

        Return:
        -------
            : Bool
            True if Parmed object deletion was successful
        """

        if not self.has_parmed:
            raise ValueError("The Parmed structure is not present on the record")

        self.rec.delete_field(Fields.pmd_structure)

        return True

    def __getattr__(self, name):
        try:
            return getattr(self.rec, name)
        except AttributeError:
            raise AttributeError(
                "'%s' object has no attribute '%s'" % (type(self).__name__, name))
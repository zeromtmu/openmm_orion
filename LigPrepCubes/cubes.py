# (C) 2018 OpenEye Scientific Software Inc. All rights reserved.
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

import traceback
from LigPrepCubes import ff_utils
from floe.api import ParallelMixin, parameter

from cuberecord import OERecordComputeCube

from cuberecord.ports import  RecordOutputPort

from oeommtools import utils as oeommutils

from Standards import Fields

from floe.constants import ADVANCED

from openeye import oechem


class LigandChargeCube(ParallelMixin, OERecordComputeCube):
    title = "Ligand Charge Cube"
    version = "0.0.0"
    classification = [["Ligand Preparation", "OEChem", "Ligand preparation"]]
    tags = ['OEChem', 'Quacpac']
    description = """
    This cube charges ligands by using the ELF10 charge method. If the ligands
    are already charged the cube parameter charge_ligand can be used to skip the
    charging stage

    Input:
    -------
    oechem.OEMCMol - Streamed-in of the ligand molecules

    Output:
    -------
    oechem.OEMCMol - Emits the charged ligands
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    max_conformers = parameter.IntegerParameter(
        'max_conformers',
        default=800,
        help_text="Max number of ligand conformers generated to charge the ligands")

    charge_ligands = parameter.BooleanParameter(
        'charge_ligands',
        default=True,
        description='Flag used to set if charge the ligands or not')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:

            if not record.has_value(Fields.primary_molecule):
                raise ValueError("Missing Primary Molecule field")

            ligand = record.get_value(Fields.primary_molecule)

            if oechem.OECalculateMolecularWeight(ligand) > 900.0:  # Units are in Dalton
                self.opt['Logger'].warn("[{}] The molecule {} seems to have a large molecular "
                                        "weight for a ligand: {:.2f} Da"
                                        .format(self.title,
                                                ligand.GetTitle(),
                                                oechem.OECalculateMolecularWeight(ligand)))

            # Ligand sanitation
            ligand = oeommutils.sanitizeOEMolecule(ligand)

            # Charge the ligand
            if self.opt['charge_ligands']:
                charged_ligand = ff_utils.assignELF10charges(ligand,
                                                             self.opt['max_conformers'],
                                                             strictStereo=False)

                # If the ligand has been charged then transfer the computed
                # charges to the starting ligand
                map_charges = {at.GetIdx(): at.GetPartialCharge() for at in charged_ligand.GetAtoms()}
                for at in ligand.GetAtoms():
                    at.SetPartialCharge(map_charges[at.GetIdx()])
                self.log.info("[{}] ELF10 charge method applied to the ligand: {}".format(self.title,
                                                                                          ligand.GetTitle()))

            record.set_value(Fields.primary_molecule, ligand)

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed record
            self.failure.emit(record)


class LigandSetting(OERecordComputeCube):
    title = "Ligand Setting"
    version = "0.0.0"
    classification = [["Ligand Preparation", "OEChem", "Ligand preparation"]]
    tags = ['OEChem']
    description = """
    This cube is setting the ligand to perform MD. The ligand residue name is set as
    LIG, each ligand is identified with an integer number and a title is generated by
    using the ligand name 

    Input:
    -------
    Data record Stream - Streamed-in of the ligand molecules

    Output:
    -------
    Data Record Stream - Emits the MD set ligands
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    # Ligand Residue Name
    lig_res_name = parameter.StringParameter('lig_res_name',
                                             default='LIG',
                                             help_text='The ligand residue name')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.count = 0

    def process(self, record, port):
        try:
            if not record.has_value(Fields.primary_molecule):
                self.log.error("Missing '{}' field".format(Fields.primary_molecule.get_name()))
                raise ValueError("Missing Primary Molecule")

            ligand = record.get_value(Fields.primary_molecule)

            if ligand.NumConfs() > 1:
                self.opt['Logger'].info("[{}] The molecule {} has multiple conformers. Each single conformer "
                                        "will be treated as a new molecule".format(self.title,
                                                                                   ligand.GetTitle()))

            if oechem.OECalculateMolecularWeight(ligand) > 900.0:  # Units are in Dalton
                self.opt['Logger'].warn("[{}] The molecule {} seems to have a large molecular "
                                        "weight for a ligand: {:.2f} Da"
                                        .format(self.title,
                                                ligand.GetTitle(),
                                                oechem.OECalculateMolecularWeight(ligand)))

            for at in ligand.GetAtoms():
                residue = oechem.OEAtomGetResidue(at)
                residue.SetName(self.args.lig_res_name)
                oechem.OEAtomSetResidue(at, residue)

            num_conf_counter = 0

            for conf in ligand.GetConfs():

                conf_mol = oechem.OEMol(conf)

                name = ligand.GetTitle()[0:12]

                if not name:
                    name = 'LIG' + str(self.count)

                ligand_title = 'l' + name

                if ligand.GetMaxConfIdx() > 1:
                    ligand_title += '_c' + str(num_conf_counter)

                # conf_mol.SetTitle(ligand_title)

                record.set_value(Fields.id, self.count)
                record.set_value(Fields.title, ligand_title)
                record.set_value(Fields.primary_molecule, conf_mol)

                num_conf_counter += 1

                self.count += 1

                self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed record
            self.failure.emit(record)
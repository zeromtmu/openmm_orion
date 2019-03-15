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

import traceback
from floe.api import (ParallelMixin,
                      parameter)

from cuberecord import OERecordComputeCube

from oeommtools import utils as oeommutils

from MDOrion.ForceField import utils as ffutils

import parmed

from openeye import oechem

from oeommtools import data_utils as pack_utils

from MDOrion.Standards import (MDStageTypes,
                               Fields,
                               MDRecords)

from MDOrion.MDEngines.utils import MDState

from simtk.openmm import app
from simtk import unit


class ForceFieldCube(ParallelMixin, OERecordComputeCube):
    title = "Force Field Application"
    version = "0.1.0"
    classification = [["Force Field"]]
    tags = ['ForceField']
    description = """
    This cube parametrized a system with the selected force fields. 
    The cube tries to split a system into components: protein, ligand, 
    water and excipients. The user can select the parametrization to be 
    applied to each component. The protein forcefield is limited to 
    standard amino acids and limited support to non-standard. Sugars 
    are not currently supported but this will be improved in coming 
    releases. The cube requires a record as input and produces a new 
    record where the system has been parametrized. The parametrization 
    is carried out by using a Parmed object 
    (https://github.com/ParmEd/ParmEd) 
    which will be present on the emitted record. The supported protein 
    force fields are amber99sb-ildn and the new amberfb-15. Small organic
    molecules like ligands and excipients can be parametrized by using 
    GAFF, GAFF2 and SMIRNOFF. The system spitting is based on the ligand 
    residue name. The default one is “LIG” and can be changed by using 
    the provided cube parameter. Water is currently parametrized by 
    using TIP3P force field water model only.

    Input:
    -------
    oechem.OEDataRecord - Streamed-in of the systems to parametrize

    Output:
    -------
    oechem.OEDataRecord - Streamed-out of records with the parametrized systems.
    Each record will contain a new Parmed object that carry out the 
    system parametrization
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    protein_forcefield = parameter.StringParameter(
        'protein_forcefield',
        default='amber99sbildn.xml',
        choices=['amber99sbildn.xml', 'amberfb15.xml'],
        help_text='Force field parameters to be applied to the protein')

    solvent_forcefield = parameter.StringParameter(
        'solvent_forcefield',
        default='tip3p.xml',
        help_text='Force field parameters to be applied to the water')

    ligand_forcefield = parameter.StringParameter(
        'ligand_forcefield',
        default='GAFF',
        choices=['GAFF', 'GAFF2', 'SMIRNOFF'],
        help_text='Force field to be applied to the ligand')

    lig_res_name = parameter.StringParameter(
        'lig_res_name',
        default='LIG',
        help_text='Ligand residue name. This is used during the spitting to identify the ligand')

    other_forcefield = parameter.StringParameter(
        'other_forcefield',
        default='GAFF',
        choices=['GAFF', 'GAFF2', 'SMIRNOFF'],
        help_text='Force field used to parametrize other molecules not recognized by the '
                  'protein force field like excipients')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            opt = self.opt
            opt['CubeTitle'] = self.title

            if not record.has_value(Fields.primary_molecule):
                self.log.error("Missing molecule Primary Molecule' field")
                self.failure.emit(record)
                return

            system = record.get_value(Fields.primary_molecule)

            if not record.has_value(Fields.title):
                self.log.warn("Missing record Title field")
                system_title = system.GetTitle()[0:12]
            else:
                system_title = record.get_value(Fields.title)

            # Split the complex in components in order to apply the FF
            protein, ligand, water, excipients = oeommutils.split(system, ligand_res_name=opt['lig_res_name'])

            self.log.info("[{}] \nComplex name: {}\nProtein atom numbers = {}\nLigand atom numbers = {}\n"
                          "Water atom numbers = {}\nExcipients atom numbers = {}".format(opt['CubeTitle'],
                                                                                         system_title,
                                                                                         protein.NumAtoms(),
                                                                                         ligand.NumAtoms(),
                                                                                         water.NumAtoms(),
                                                                                         excipients.NumAtoms()))
            if not record.has_value(Fields.id):
                raise ValueError("Missing ID Field")

            sys_id = record.get_value(Fields.id)

            # Unique prefix name used to output parametrization files
            opt['prefix_name'] = system_title + '_'+str(sys_id)

            oe_mol_list = []
            par_mol_list = []

            # Apply FF to the Protein
            if protein.NumAtoms():
                oe_mol_list.append(protein)
                protein_structure = ffutils.applyffProtein(protein, opt)
                par_mol_list.append(protein_structure)

            # Apply FF to the ligand
            if ligand.NumAtoms():
                oe_mol_list.append(ligand)
                ligand_structure = ffutils.applyffLigand(ligand, opt)
                par_mol_list.append(ligand_structure)

            # Apply FF to water molecules
            if water.NumAtoms():
                oe_mol_list.append(water)
                water_structure = ffutils.applyffWater(water, opt)
                par_mol_list.append(water_structure)

            # Apply FF to the excipients
            if excipients.NumAtoms():
                excipient_structure = ffutils.applyffExcipients(excipients, opt)
                par_mol_list.append(excipient_structure)

                # The excipient order is set equal to the order in related
                # parmed structure to avoid possible atom index mismatching
                excipients = oeommutils.openmmTop_to_oemol(excipient_structure.topology,
                                                           excipient_structure.positions,
                                                           verbose=False)
                oechem.OEPerceiveBondOrders(excipients)
                oe_mol_list.append(excipients)

            # Build the overall Parmed structure
            system_structure = parmed.Structure()

            for struc in par_mol_list:
                system_structure = system_structure + struc

            system_reassembled = oe_mol_list[0].CreateCopy()
            num_atom_system = system_reassembled.NumAtoms()

            for idx in range(1, len(oe_mol_list)):
                oechem.OEAddMols(system_reassembled, oe_mol_list[idx])
                num_atom_system += oe_mol_list[idx].NumAtoms()

            if not num_atom_system == system_structure.topology.getNumAtoms():
                raise ValueError("Parmed and OE topologies mismatch atom number "
                                 "error for system: {}".format(system_title))

            system_reassembled.SetTitle(system.GetTitle())

            # Set Parmed structure box_vectors
            is_periodic = True

            try:
                vec_data = pack_utils.getData(system_reassembled, tag='box_vectors')
                vec = pack_utils.decodePyObj(vec_data)
                system_structure.box_vectors = vec
            except:
                is_periodic = False
                self.log.warn("System {} has been parametrize without periodic box vectors ".format(system_title))

            # Set atom serial numbers, Ligand name and HETATM flag
            for at in system_reassembled.GetAtoms():
                thisRes = oechem.OEAtomGetResidue(at)
                thisRes.SetSerialNumber(at.GetIdx())
                if thisRes.GetName() == 'UNL':
                    # thisRes.SetName("LIG")
                    thisRes.SetHetAtom(True)
                oechem.OEAtomSetResidue(at, thisRes)

            if system_reassembled.GetMaxAtomIdx() != system_structure.topology.getNumAtoms():
                raise ValueError("OEMol system {} and generated Parmed structure "
                                 "mismatch atom numbers".format(system_title))

            # Copying the charges between the parmed structure and the oemol
            for parm_at, oe_at in zip(system_structure.atoms, system_reassembled.GetAtoms()):

                if parm_at.atomic_number != oe_at.GetAtomicNum():
                    raise ValueError("Atomic number mismatch between the Parmed and the OpenEye topologies: {} - {}".
                                     format(parm_at.atomic_number, oe_at.GetAtomicNum()))

                oe_at.SetPartialCharge(parm_at.charge)

            # Check if it is possible to create the OpenMM System
            if is_periodic:
                omm_system = system_structure.createSystem(nonbondedMethod=app.CutoffPeriodic,
                                                           nonbondedCutoff=10.0 * unit.angstroms,
                                                           constraints=None,
                                                           removeCMMotion=False,
                                                           rigidWater=False)
            else:
                omm_system = system_structure.createSystem(nonbondedMethod=app.NoCutoff,
                                                           constraints=None,
                                                           removeCMMotion=False,
                                                           rigidWater=False)

            record.set_value(Fields.primary_molecule, system_reassembled)
            record.set_value(Fields.title, system_title)

            record.set_value(Fields.pmd_structure, system_structure)

            mdstate = MDState(system_structure)

            md_stage = MDRecords.MDStageRecord(self.title, MDStageTypes.SETUP,
                                               MDRecords.MDSystemRecord(system_reassembled, mdstate))

            record.set_value(Fields.md_stages, [md_stage])

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed record
            self.failure.emit(record)

        return

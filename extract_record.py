from Standards import Fields
from cuberecord.converters.oldrecordutil import oe_mol_to_data_record
from openeye import oechem
from datarecord import OEField, Types

mol_stream = oechem.oemolistream("prep.oeb")
mol = oechem.OEMol()

oechem.OEReadMolecule(mol_stream, mol)
mol_stream.close()

record = oe_mol_to_data_record(mol)

stages = record.get_value(Fields.md_stages)
print(stages)
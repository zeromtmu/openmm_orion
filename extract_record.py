from Standards import Fields
from cuberecord.converters.oldrecordutil import oe_mol_to_data_record
from openeye import oechem
from datarecord import OEField, Types
from cuberecord.cube_testing import OEMolRecordStream


# mol_stream = oechem.oemolistream("test.oeb")
# mol = oechem.OEMol()
#
# oechem.OEReadMolecule(mol_stream, mol)
# mol_stream.close()
#
# record = oe_mol_to_data_record(mol)
#
# stages = record.get_value(Fields.md_stages)
#
# stage = stages[-1]
#
# print(stage)


ifs = OEMolRecordStream("test.oedb")

for record in ifs:
    stages = record.get_value(Fields.md_stages)
    stage = stages[-1]
    print(stage.get_value(Fields.log_data))

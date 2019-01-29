
from openeye import oechem

from datarecord import read_mol_record

from Standards import Fields

filename = "prep.oedb"

ifs = oechem.oeifstream(filename)

records = []

while True:
    record = read_mol_record(ifs)
    if record is None:
        break
    records.append(record)
ifs.close()

for record in records:
    pmd = record.get_value(Fields.pmd_structure)
    pmd.save("system.pdb")


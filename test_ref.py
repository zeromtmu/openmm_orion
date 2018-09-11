from openeye import oechem

from Standards import Fields

from datarecord import read_mol_record

ifs = oechem.oeifstream("prep3.oedb")
records = []
while True:
    record = read_mol_record(ifs)
    if record is None:
        break
    records.append(record)
ifs.close()


for record in records:
    #pmd = record.get_value(Fields.structure)
    stages = record.get_value(Fields.md_stages)
    print(len(stages))
    #stage = stages[-1]
    #mdsystem = stage.get_value(Fields.md_system)


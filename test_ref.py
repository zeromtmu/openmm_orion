from openeye import oechem

from Standards import Fields

from datarecord import read_mol_record

from MDCubes.mdutils import MDStructure

ifs = oechem.oeifstream("prep.oedb")
records = []
while True:
    record = read_mol_record(ifs)
    if record is None:
        break
    records.append(record)
ifs.close()


for record in records:
    stages = record.get_value(Fields.md_stages)
    stage = stages[-1]
    mdsystem = stage.get_value(Fields.md_system)
    pmd = mdsystem.get_value(Fields.structure)


mdstruct = MDStructure(pmd)

state =mdstruct.get_state()

print(state.velocities)
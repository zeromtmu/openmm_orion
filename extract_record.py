from Standards import Fields
from tempfile import NamedTemporaryFile
from openeye import oechem
from datarecord import OEField, Types
from cuberecord.cube_testing import OEMolRecordStream
from datarecord import OEReadRecords
from datarecord import OEWriteRecord

from orionclient.types import Dataset
from orionclient.session import APISession

from datarecord import read_mol_record

from datarecord import OEWriteRecord

from simtk.openmm import app

# ds = Dataset.upload(APISession, "foobar", "p38_l38_a_2n_nvt_5ns.oeb.gz")
from oeommtools import utils as oeommutils

# ifs = OEMolRecordStream("p38_l38_a_2n_nvt_5ns.oeb.gz")

ifs = oechem.oeifstream("pP38_lig38a_2n_npt_5ns.oedb")
records = []
while True:
    record = read_mol_record(ifs)
    if record is None:
        break
    records.append(record)
ifs.close()

print(len(records))


for record in records:
    # stages = record.get_value(Fields.md_stages)
    # print("Len stages = {}".format(len(stages)))
    # stage = stages[-1]
    # #print(stage.has_value(Fields.log_data))
    # #print(stage.get_value(Fields.log_data))
    # mdsystem = stage.get_value(Fields.md_system)
    pmd = record.get_value(Fields.pmd_structure)
    pmd_split = pmd.split()[2][0]
    pmd_split.save("cazzo.pdb")
    print(pmd_split)
    #pmd.save("sys.pdb", overwrite=True)

    # complex = record.get_value(Fields.primary_molecule)
    #
    # # Split the complex in components in order to apply the FF
    # protein, ligand, water, excipients = oeommutils.split(complex, ligand_res_name='LIG')
    #
    # print("Protein atom numbers = {}\nLigand atom numbers = {}\n"
    #               "Water atom numbers = {}\nExcipients atom numbers = {}".format(
    #                                                                              protein.NumAtoms(),
    #                                                                              ligand.NumAtoms(),
    #                                                                              water.NumAtoms(),
    #                                                                              excipients.NumAtoms()))
    #
    # topology, positions = oeommutils.oemol_to_openmmTop(water)
    #
    # ff = app.ForceField('amber99sbildn.xml', 'tip4pew.xml')
    #
    # modeller = app.Modeller(topology, positions)
    # modeller.addExtraParticles(ff)
    #
    # app.PDBFile.writeFile(modeller.topology, modeller.positions, open('tip4pew.pdb', 'w'))


    # complex = mdsystem.get_value(Fields.topology)
    # with oechem.oemolostream("compl.oeb") as ofs:
    #     oechem.OEWriteConstMolecule(ofs, complex)
    # ligand = record.get_value(Fields.ligand)
    # with oechem.oemolostream("ligand.oeb") as ofs:
    #     oechem.OEWriteConstMolecule(ofs, ligand)

    #
# ofs = oechem.oeofstream("pP38_lig38a_2n_nvt_5ns_mod.oedb")
#
# for record in records:
#     stages = record.get_value(Fields.md_stages)
#     stage = stages[-1]
#     stages = [stage]
#     record.set_value(Fields.md_stages, stages)
#     OEWriteRecord(ofs, record, fmt='binary')
#
# ofs.close()

#for record in ifs:
    #print([(x.get_name(), x.get_type()) for x in record.get_fields()])
    # print(record.get_value(Fields.title))
    # for field in record.get_fields():
    #     print(field.get_name())
#     #print(record.get_value(OEField("DG", Types.Float)))
    #assert record.has_value(Fields.md_stages)
        # raise ValueError("The System does not seem to be parametrized by the Force Field")
    # stages = record.get_value(Fields.md_stages)
    # print("Len stages = {}".format(len(stages)))
    # stage = stages[-1]
    # #print(stage.get_value(Fields.trajectory))
    # print(stage.has_value(Fields.log_data))
    # print(stage.get_value(Fields.log_data))
    #print(stage.has_value(Fields.trajectory))
    #print(stage.get_value(Fields.trajectory))
    # mdsystem = stage.get_value(OEField("MDSystem", Types.Record))
    # complex = mdsystem.get_value(OEField('Topology_OEMol', Types.Chem.Mol))
    # ofs = oechem.oemolostream("complex.oeb")
    # oechem.OEWriteConstMolecule(ofs, complex)

    # for field in mdsystem.get_fields():
    #     print(field.get_name())





# for record in ifs:
#     # stages = record.get_value(Fields.md_stages)
#     # stage = stages[0]
#     # print(stage.get_value(OEField("Log_data", Types.String)))
#
#     #print(stage.get_value(Fields.log_data))
#     # if record.has_value(OEField("Title_PLMD", Types.String)):
#     #     print(record.get_value(OEField("Title_PLMD", Types.String)))
#     for field in record.get_fields():
#         print(field.get_name())
#         type_field = field.get_type()
#         print(type_field)
#     # #print(record.get_value(OEField("DG", Types.Float)))
#     # protein = record.get_value(OEField("Protein", Types.Chem.Mol))
#     # ofs = oechem.oemolostream("protein"+str(count)+".oeb")
#     # oechem.OEWriteConstMolecule(ofs, protein)
#     count += 1

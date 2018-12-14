
from openeye import oechem

from openeye import oedepict

ifs = oechem.oemolistream("freesolv_mini.oeb.gz")

for mol in ifs.GetOEGraphMols():
    oedepict.OEPrepareDepiction(mol)
    width, height = 300, 300
    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    disp = oedepict.OE2DMolDisplay(mol, opts)
    oedepict.OERenderMolecule(mol.GetTitle()+'.svg', disp)

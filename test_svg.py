
from openeye import oechem

from openeye import oedepict

ifs = oechem.oemolistream("freesolv_mini.oeb")

for mol in ifs.GetOEGraphMols():
    oedepict.OEPrepareDepiction(mol)
    oedepict.OERenderMolecule(mol.GetTitle()+'.svg', mol)

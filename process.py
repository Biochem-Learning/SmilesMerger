from rdkit import Chem
from rdkit.Chem import rdmolops

def combine_smile(smile1, smile2, idx1=0, idx2=0):
    mol1 = Chem.MolFromSmiles(smile1)
    mol2 = Chem.MolFromSmiles(smile2)

    # Combine the two molecules
    combo = rdmolops.CombineMols(mol1, mol2)

    # Convert to editable mol
    editable = Chem.RWMol(combo)

    # mol1 has N atoms; mol2's atoms are offset by N
    offset = mol1.GetNumAtoms()

    # Add a bond between mol1[idx1] and mol2[idx2]
    editable.AddBond(idx1, offset + idx2, Chem.rdchem.BondType.SINGLE)

    # Final molecule
    final = editable.GetMol()

    # Sanitize and convert to SMILES
    Chem.SanitizeMol(final)
    return Chem.MolToSmiles(final)

    # convert smile string to mol
    # combine 2 mol
    # convert back to smilestr
    # return smilestr
    # get number of explicit H
    # set number of explicat H
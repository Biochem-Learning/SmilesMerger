from rdkit import Chem
from rdkit.Chem import rdmolops

def combine_smile(smile1, smile2):
    mol1 = Chem.MolFromSmiles(smile1)
    mol2 = Chem.MolFromSmiles(smile2)

    if mol1 is None or mol2 is None:
        raise ValueError("One or both SMILES strings are invalid.")

    combined_mol = Chem.ReplaceSubstructs(mol1, Chem.MolFromSmarts("C"), mol2)
    combined_smiles = Chem.MolToSmiles(combined_mol[0])

    print(combined_smiles)
    return combined_smiles

    # convert smile string to mol
    # combine 2 mol
    # convert back to smilestr
    # return smilestr
    # get number of explicit H
    # set number of explicat H

combine_smile("CC(C)=O", "[H]")
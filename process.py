from rdkit import Chem
from rdkit.Chem import rdmolops

def combine_smile(smile1, smile2, connection_point_index=1):
    mol1 = Chem.MolFromSmiles(smile1)
    mol2 = Chem.MolFromSmiles(smile2)
    query_smiles = "C"
    query = Chem.MolFromSmiles(query_smiles) 

    if mol1 is None or mol2 is None:
        raise ValueError("One or both SMILES strings are invalid.")

    placeholder_num = 0
    combined_smiles = smile1
    combined_mol = mol1

    while not is_string_at_index(combined_smiles, connection_point_index, smile2):
        combined_mol = Chem.MolFromSmiles(combined_smiles)
        combined_mol = Chem.ReplaceSubstructs(combined_mol, query, mol2)
        combined_smiles = Chem.MolToSmiles(combined_mol[0])
        placeholder_num += 1

    combined_smiles = replace_last_n(combined_smiles,  smile2, query_smiles, placeholder_num-1)
    print(combined_smiles)
    return combined_smiles

def get_reversed_index(total_length, forward_index):
    reverse_index = total_length - forward_index + 1
    return reverse_index

def replace_last_n(s, old, new, n):
    return new.join(s.rsplit(old, n))

    # convert smile string to mol
    # combine 2 mol
    # convert back to smilestr
    # return smilestr
    # get number of explicit H
    # set number of explicat H

def is_string_at_index(string1, index, string2):
    """
    Check if string2 occurs at the given index in string1.
    
    Args:
        string1 (str): The main string to search in.
        index (int): The index at which to check for string2.
        string2 (str): The substring to match.

    Returns:
        bool: True if string2 starts at index in string1, False otherwise.
    """
    return string1[index:index+len(string2)] == string2

# def combine_smile(smile1, smile2):
#     mol1 = Chem.MolFromSmiles(smile1)
#     mol2 = Chem.MolFromSmiles(smile2)

#     if mol1 is None or mol2 is None:
#         raise ValueError("One or both SMILES strings are invalid.")

#     combined_mol = Chem.ReplaceSubstructs(mol1, Chem.MolFromSmarts("C"), mol2)
#     combined_smiles = Chem.MolToSmiles(combined_mol[0])

#     print(combined_smiles)
#     return combined_smiles

combine_smile("CN(C)C", "Cl", 3)
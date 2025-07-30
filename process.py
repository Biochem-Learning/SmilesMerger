from rdkit import Chem
from rdkit.Chem import Draw

 
def mol_with_atom_index( mol ):
    atoms = mol.GetNumAtoms()
    for idx in range( atoms ):
        mol.GetAtomWithIdx( idx ).SetProp( 'molAtomMapNumber', str( mol.GetAtomWithIdx( idx ).GetIdx() ) )
    return mol

# def combine_smile(smile1, smile2, connection_point_index=1):
#     mol1 = Chem.MolFromSmiles(smile1)
#     mol2 = Chem.MolFromSmiles(smile2)

#     if mol1 is None or mol2 is None:
#         raise ValueError("One or both SMILES strings are invalid.")

#     combined_mol = Chem.CombineMols(mol1,mol2)

#     final_mol = mol_with_atom_index(combined_mol)

#     ed_combined_mol = Chem.EditableMol(combined_mol)

#     ed_combined_mol.AddBond(1, 4, order=Chem.rdchem.BondType.SINGLE)
#     # Draw with atom numbers

#     final_mol1 = ed_combined_mol.GetMol()
#     return Draw.MolToImage(final_mol1, size=(400, 300))


# def combine_smile(smile1, smile2, connection_point_index=1):
#     mol1 = Chem.MolFromSmiles(smile1)
#     mol2 = Chem.MolFromSmiles(smile2)
#     query_smiles = "C"
#     query = Chem.MolFromSmiles(query_smiles) 

#     if mol1 is None or mol2 is None:
#         raise ValueError("One or both SMILES strings are invalid.")

#     placeholder_num = 0
#     combined_smiles = smile1
#     combined_mol = mol1

#     while not is_string_at_index(combined_smiles, connection_point_index, smile2):
#         combined_mol = Chem.MolFromSmiles(combined_smiles)
#         combined_mol = Chem.ReplaceSubstructs(combined_mol, query, mol2)
#         combined_smiles = Chem.MolToSmiles(combined_mol[0])
#         placeholder_num += 1

#     combined_smiles = replace_last_n(combined_smiles,  smile2, query_smiles, placeholder_num-1)
#     print(combined_smiles)
#     return combined_smiles

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

#     combined_mol = Chem.ReplaceSubstructs(mol1, Chem.MolFromSmiles("C"), mol2)[0]
#     combined_smiles = Chem.MolToSmiles(combined_mol)

#     print(combined_smiles)
#     return Draw.MolToImage(combined_mol)

def combine_smile(smile1, smile2, connection_point_index=1):
    mol1 = Chem.MolFromSmiles(smile1)
    mol2 = Chem.MolFromSmiles(smile2)

    if mol1 is None or mol2 is None:
        raise ValueError("One or both SMILES strings are invalid.")

    # index_mol1 = mol_with_atom_index(mol1)

    combined_mol = Chem.ReplaceSubstructs(mol1, Chem.MolFromSmarts("[CH3]"), mol2)[0]
    combined_smiles = Chem.MolToSmiles(combined_mol)
    
    # reconstructed_smiles = reconstruct_smiles_based_on_atom_idx(combined_smiles, smile2)

    print(combined_smiles)
    return combined_smiles

    # reconstruct combined_smiles

    # replace all the element in smile with placeholder
    # based on index, insert atom into smile
import re

def reconstruct_smiles_based_on_atom_idx(smiles,frag2):
    placeholdered_smiles = replace_atoms_with_placeholder(smiles)

    atom_dict = dict_of_indexed_atom(smiles)

    partial_reconstructed_smiles= replace_stars_with_atoms(placeholdered_smiles, atom_dict)
    reconstructed_smiles = partial_reconstructed_smiles.replace("*", frag2)
    return reconstructed_smiles

def replace_stars_with_atoms(placeholdered_smiles, atom_dict):
    result = ""
    index = 0
    i = 0

    while i < len(placeholdered_smiles):
        if placeholdered_smiles[i] == "*":
            result += atom_dict.get(index)  # fallback to "*" if index missing
            index += 1
            i += 1
        else:
            result += placeholdered_smiles[i]
            i += 1
    return result

def replace_atoms_with_placeholder(smiles):
    result = ""
    i = 0
    while i < len(smiles):
        if smiles[i] == "[":
            # Find the closing bracket
            j = i
            while j < len(smiles) and smiles[j] != "]":
                j += 1
            if j < len(smiles):
                result += "*"
                i = j + 1
            else:
                result += smiles[i]
                i += 1
        elif smiles[i].isalpha():
            # If it's a two-letter element like Cl, Br, etc.
            if i + 1 < len(smiles) and smiles[i + 1].islower():
                result += "*"
                i += 2
            else:
                result += "*"
                i += 1
        else:
            result += smiles[i]
            i += 1
    return result

def dict_of_indexed_atom(smiles):
    mol = Chem.MolFromSmiles(smiles)
    full_smiles = smiles

    atom_fragments = {}
    seen_fragments = set()

    for i in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(i)
        atom_map_num = atom.GetAtomMapNum()

        match = re.search(rf"\[[^\]]+:{atom_map_num}\]", full_smiles)
        if match:
            atom_fragment = match.group(0)
        else:
            atom_fragment = atom.GetSymbol()

        # Avoid duplicate fragments
        if atom_fragment in seen_fragments:
            atom_fragment = "*"
        else:
            seen_fragments.add(atom_fragment)

        atom_fragments[i] = atom_fragment

    return atom_fragments 

# def attach_frag_to_atom(mol1_smiles, mol2_smiles, mol1_attach_idx):
#     mol1 = Chem.MolFromSmiles(mol1_smiles)
#     mol2 = Chem.MolFromSmiles(mol2_smiles)

#     # 1. Find the atom in mol2 that has only 1 neighbor — this will be the connection point
#     if mol2.GetNumAtoms() == 1:
#         mol2_attach_idx = 0
#     else:
#         mol2_attach_idx = None
#         for atom in mol2.GetAtoms():
#             if len(atom.GetNeighbors()) == 1:
#                 mol2_attach_idx = atom.GetIdx()
#                 break
#         if mol2_attach_idx is None:
#             raise ValueError("No terminal atom found in mol2.")

#     # 2. Add dummy atoms to use as placeholders for bonding
#     rw_mol1 = Chem.RWMol(mol1)
#     rw_mol2 = Chem.RWMol(mol2)

#     rw_mol1.AddAtom(Chem.Atom(0))  # Dummy atom with atomic number 0
#     dummy_idx1 = rw_mol1.GetNumAtoms() - 1
#     rw_mol1.AddBond(mol1_attach_idx, dummy_idx1, Chem.BondType.SINGLE)

#     rw_mol2.AddAtom(Chem.Atom(0))
#     dummy_idx2 = rw_mol2.GetNumAtoms() - 1
#     rw_mol2.AddBond(mol2_attach_idx, dummy_idx2, Chem.BondType.SINGLE)

#     # 3. Combine molecules via CombineMols and use editable mol
#     combo = Chem.CombineMols(rw_mol1, rw_mol2)
#     edcombo = Chem.EditableMol(Chem.RWMol(combo))

#     # Calculate offset
#     offset = rw_mol1.GetNumAtoms()

#     # 4. Remove dummy atoms and add bond
#     edcombo.AddBond(mol1_attach_idx, mol2_attach_idx + offset, Chem.BondType.SINGLE)
#     edcombo.RemoveAtom(dummy_idx2 + offset)
#     edcombo.RemoveAtom(dummy_idx1)
#     merged = edcombo.GetMol()

#     # 5. Sanitize and return SMILES
#     Chem.SanitizeMol(merged)
#     return Chem.MolToSmiles(merged, canonical=False)

# mol1_smiles = "CC(C)(C)Cl"
# mol2_smiles = "Cl"

# # combine_smile(mol1_smiles, mol2_smiles, 4)
# print(attach_frag_to_atom("CC(C)(C)Cl", "CN(C)C", 3))


# CC(C)(C)C
# COC
# [H]C
# ClC
# BrC
# CC(C)=O
# CN(C)C
# CC(C)=C(C)C
# CC(=C(C)C(C)=C1C)C(C)=C1C
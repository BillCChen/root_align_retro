# split = "train"
# split_num = 0
# base_dir = f"/root/reaction_data/pretrain_aug/CHORISO_PtoTMPtoR_aug5/{split}/raw_mapping_reaction-{split}-merged.txt"
# save_dir = f"/root/reaction_data/pretrain_aug/CHORISO_PtoTMPtoR_aug5/{split}/tmp_smarts_{split_num}.txt"
# error_dir = f"/root/reaction_data/pretrain_aug/CHORISO_PtoTMPtoR_aug5/{split}/error_reaction_{split_num}.txt"
# import numpy as np
# import pandas as pd
# import argparse
# import os
# import re
# import random
# import textdistance
# import multiprocessing

# from rdkit import Chem
# from tqdm import tqdm


# from rdkit import RDLogger
# RDLogger.DisableLog('rdApp.*')
# import sys
# sys.path.append("/root/retro_synthesis/template_analysis")
# # 添加 /root/retro_synthesis/reaction_utils 的库查找路径
# sys.path.append("/root/retro_synthesis/reaction_utils")
# from rxnutils.chem.reaction import ChemicalReaction
# def smi_tokenizer(smi):
#     pattern = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
#     regex = re.compile(pattern)
#     tokens = [token for token in regex.findall(smi)]
#     assert smi == ''.join(tokens)
#     return ' '.join(tokens)

# def precise_tempalte_extraction(reaction_smiles,jump=3,expand_ring=True):
#     template =  ChemicalReaction(reaction_smiles).generate_reaction_template(jump,expand_ring)[0].smarts
#     # 确保不要浪费内存
#     assert type(template) == str, "Template extraction failed, expected a string."
#     return template


# def clear_map_canonical_smiles(smi, canonical=True, root=-1, type_='smarts'):
#     assert type_ in ['smiles', 'smarts'], "Type must be either 'smiles' or 'smarts'."
#     if type_ == 'smiles':
#         mol = Chem.MolFromSmiles(smi,sanitize=False)
#     elif type_ == 'smarts':
#         mol = Chem.MolFromSmarts(smi)
#     else:
#         raise ValueError("Unsupported type: {}".format(type_))
#     if mol is not None:
#         for atom in mol.GetAtoms():
#             if atom.HasProp('molAtomMapNumber'):
#                 atom.ClearProp('molAtomMapNumber')
#         return Chem.MolToSmiles(mol, isomericSmiles=True, rootedAtAtom=root, canonical=canonical)
#     else:
#         return smi
# def remove_atom_mapping_template(template):
#     """
#     Remove atom mapping from a reaction template.
#     Args:
#         template (str): Reaction template with atom mapping.
#         mode: 
#             single: only use function_clear_map_canonical_smiles
#             double: use function_clear_map_canonical_smiles and Chem.MolToSmiles(Chem.MolFromSmiles()
#     Returns:
#         str: Reaction template without atom mapping.
#     """
#     def remove_mapping(smi):
#         return clear_map_canonical_smiles(smi, canonical=True, type_='smarts')

#     left , right = template.split('>>')
#     if len(left.split('.')) > 1:
#         left = '.'.join([remove_mapping(c) for c in left.split('.')])
#     else:
#         left = remove_mapping(left)
#     if len(right.split('.')) > 1:
#         right = '.'.join([remove_mapping(c) for c in right.split('.')])
#     else:
#         right = remove_mapping(right)
#     return left + '>>' + right
# with open(base_dir, 'r') as f:
#     data = f.readlines()
# data = [line.strip() for line in data if line.strip()]
# # print(f"Processing {len(data)} from {split_num*400000} to {(split_num+1)*400000} reactions.")
# # data = data[split_num*400000 : split_num*400000 + 400000]

# all_tmp = []
# store_temp = ""
# from tqdm import tqdm
# error_num = 0
# error_reaction = []
# for reaction_ in tqdm(data, desc="Processing reactions"):
# # for reaction_ in data:
#     reactant , product, num = reaction_.split('>>')

#     if num == "0":
#         # print(f"{reactant}>>{product}")
#         try:
#             store_temp = precise_tempalte_extraction(f"{reactant}>>{product}",3)
#         except Exception as e:
#             print(f"Error processing reaction {reaction_[:66]}")
#             error_reaction.append(reaction_)
#             error_num += 1
#             if error_num % 100 == 0:
#                 print(f"Processed {len(all_tmp)} reactions, encountered {error_num} errors so far.")
#             continue
#         tmp_reactant,tmp_prodcut = store_temp.split('>>')
#         assert len(tmp_reactant) > 0, "The reactant part of the template cannot be empty."
#         assert len(tmp_prodcut) > 0, "The product part of the template cannot be empty."
#         all_tmp.append(store_temp)
# # assert len(all_tmp) * 5 == len(data), "The number of reactants and products must match the number of reactions."
# # 输出error_num 的数字以及总数据的占比
# print(f"Total reactions processed: {len(data)}")
# print(f"Total templates extracted: {len(all_tmp)}")
# print(f"Total errors encountered: {error_num}")
# # 输出错误占比
# print(f"Error rate: {error_num / len(data) * 100:.2f}%")
# with open(save_dir, "w") as f:
#     for smi in all_tmp:
#         f.write(smi + '\n')
# with open(error_dir, "w") as f:
#     for smi in error_reaction:
#         f.write(smi + '\n')


# # This script is for multi-processing the extraction of reaction templates from a dataset of reactions.



import numpy as np
import pandas as pd
import argparse
import os
import re
import random
import textdistance
import multiprocessing
from rdkit import Chem
from tqdm import tqdm
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
import sys
from rxnutils.chem.reaction import ChemicalReaction

def smi_tokenizer(smi):
    pattern = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    assert smi == ''.join(tokens)
    return ' '.join(tokens)

def precise_tempalte_extraction(reaction_smiles, jump=3, expand_ring=True):
    template = ChemicalReaction(reaction_smiles).generate_reaction_template(jump, expand_ring=True)[0].smarts
    assert type(template) == str, "Template extraction failed, expected a string."
    return template

def clear_map_canonical_smiles(smi, canonical=True, root=-1, type_='smarts'):
    assert type_ in ['smiles', 'smarts'], "Type must be either 'smiles' or 'smarts'."
    if type_ == 'smiles':
        mol = Chem.MolFromSmiles(smi, sanitize=False)
    elif type_ == 'smarts':
        mol = Chem.MolFromSmarts(smi)
    else:
        raise ValueError("Unsupported type: {}".format(type_))
    if mol is not None:
        for atom in mol.GetAtoms():
            if atom.HasProp('molAtomMapNumber'):
                atom.ClearProp('molAtomMapNumber')
        return Chem.MolToSmiles(mol, isomericSmiles=True, rootedAtAtom=root, canonical=canonical)
    else:
        return smi

def remove_atom_mapping_template(template):
    def remove_mapping(smi):
        return clear_map_canonical_smiles(smi, canonical=True, type_='smarts')
    left, right = template.split('>>')
    if len(left.split('.')) > 1:
        left = '.'.join([remove_mapping(c) for c in left.split('.')])
    else:
        left = remove_mapping(left)
    if len(right.split('.')) > 1:
        right = '.'.join([remove_mapping(c) for c in right.split('.')])
    else:
        right = remove_mapping(right)
    return left + '>>' + right
jump_num = 3
def process_reaction(tuple_data):
    reaction_, cls = tuple_data
    try:
        reactant, product, num = reaction_.split('>>')
        if num == "0":
            template = precise_tempalte_extraction(f"{reactant}>>{product}", jump_num, expand_ring=True)
            tmp_reactant, tmp_product = template.split('>>')
            assert len(tmp_reactant) > 0, "The reactant part of the template cannot be empty."
            assert len(tmp_product) > 0, "The product part of the template cannot be empty."
            return template, None , cls
        return None, reaction_ , None
    except Exception as e:
        return None, reaction_ , None

def main():
    # conda activate unirxn_G
    split = "train"
    split_num = ""
    base_dir = f"/root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/{split}/raw_mapping_reaction-{split}.txt"
    cls_dir = f"/root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/{split}/reaction_cls-{split}.txt"
    save_dir = f"/root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/{split}/tmp_smarts_{split_num}_jump{jump_num}_enhance_ring.txt"
    save_cls_dir = f"/root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/{split}/tmp_smarts_cls_{split_num}_jump{jump_num}.txt"
    error_dir = f"/root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/{split}/error_reaction_{split_num}_jump{jump_num}.txt"

    # Read input data
    with open(base_dir, 'r') as f:
        data = [line.strip() for line in f.readlines() if line.strip()]
    with open(cls_dir, 'r') as f:
        cls_data = [line.strip() for line in f.readlines() if line.strip()]
    assert len(data) == len(cls_data), "The number of reactions must match the number of classes."
    # Initialize multiprocessing pool with 40 cores
    pool = multiprocessing.Pool(processes=40)
    tuple_data = list(zip(data, cls_data))
    # Process reactions in parallel with progress bar
    results = []
    for result in tqdm(pool.imap_unordered(process_reaction, tuple_data), total=len(tuple_data), desc="Processing reactions"):
        results.append(result)
    
    pool.close()
    pool.join()

    # Separate successful templates and errors
    all_tmp = [r[0] for r in results if r[0] is not None]
    all_cls = [r[2] for r in results if r[0] is not None]
    error_reaction = [r[1] for r in results if r[1] is not None]
    error_num = len(error_reaction)

    # Print statistics
    print(f"Total reactions processed: {len(data)}")
    print(f"Total templates extracted: {len(all_tmp)}")
    print(f"Total errors encountered: {error_num}")
    print(f"Error rate: {error_num / len(data) * 100:.2f}%")

    # Save results
    with open(save_dir, "w") as f:
        for smi in all_tmp:
            f.write(smi + '\n')
    with open(save_cls_dir, "w") as f:
        for smi in all_cls:
            f.write(smi + '\n')
    with open(error_dir, "w") as f:
        for smi in error_reaction:
            f.write(smi + '\n')

if __name__ == "__main__":
    main()
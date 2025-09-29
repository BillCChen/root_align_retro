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





from typing import List, Tuple, Optional, Dict, Any

from rdkit import Chem
from rdkit.Chem import AllChem
from rdchiral import template_extractor as extractor
class ReactionException(Exception):
    """Custom exception raised when failing operations on a chemical reaction"""

# def _generate_rdchiral_template(
#     reactants: List[Chem.Mol],
#     products: List[Chem.Mol],
#     radius: int = 1,
#     extra_atoms: Dict[str, Any] = None,
#     expand_hetero: bool = False,
#     ) -> str:
#         """Generate a reaction template with rdChiral.

#         :param reactants: List of reactant molecules
#         :param products: List of product molecules
#         :param radius: Template radius, defaults to 1
#         :raises ReactionException: Template generation failed: could not obtained changed atoms
#         :raises ReactionException: Template generation failed: no atoms changes
#         :raises ReactionException: Template generation failed: Timed out
#         :return: Canonical reaction template
#         """
#         changed_atoms, changed_atom_tags, err = extractor.get_changed_atoms(
#             reactants=reactants, products=products
#         )
#         # 新增逻辑：如果启用了 expand_ring，则检查并扩展芳香环
#         expand_ring = True
#         if expand_ring:
#             # 获取所有反应中心原子的原子映射编号，而不是RDKit对象
#             changed_atom_map_nums = changed_atom_tags

#             # 遍历所有反应物和产物分子
#             all_mols = reactants + products
            
#             # 遍历所有环
#             for mol in all_mols:
#                 if not mol:
#                     continue
                
#                 # 使用一个集合来存储需要扩展的环的原子映射编号，避免重复
#                 ring_atoms_to_expand = set()
                
#                 for ring_atom_indices in mol.GetRingInfo().AtomRings():
#                     # 检查该环中是否有原子与反应中心相关
#                     ring_should_be_expanded = False
                    
#                     # 先找到环中所有有原子映射编号的原子
#                     ring_map_nums = [mol.GetAtomWithIdx(idx).GetAtomMapNum() for idx in ring_atom_indices if mol.GetAtomWithIdx(idx).GetAtomMapNum()]
                    
#                     # 如果环是芳香环且有原子映射编号
#                     if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atom_indices) and ring_map_nums:
                        
#                         # 遍历环中的每个原子，检查它是否接近反应中心
#                         for ring_atom_idx in ring_atom_indices:
#                             ring_atom = mol.GetAtomWithIdx(ring_atom_idx)
                            
#                             # 检查是否与任何一个反应中心原子有短路径连接
#                             for changed_atom_map_num in changed_atom_map_nums:
#                                 # 找到对应映射编号的原子
#                                 try:
#                                     # 遍历当前分子的所有原子来找到匹配的映射编号
#                                     changed_atom = None
#                                     for atom in mol.GetAtoms():
#                                         if atom.GetAtomMapNum() == int(changed_atom_map_num):
#                                             changed_atom = atom
#                                             break
                                    
#                                     if changed_atom:
#                                         # 计算最短路径
#                                         path = Chem.rdmolops.GetShortestPath(mol, ring_atom_idx, changed_atom.GetIdx())
#                                         path_len = len(path) - 1 # 路径长度是原子数减1
                                        
#                                         # 如果路径长度在设定的半径内，则标记该环需要扩展
#                                         if path_len <= radius:
#                                             ring_should_be_expanded = True
#                                             break
#                                 except Exception:
#                                     # 如果路径计算失败，跳过
#                                     continue
                            
#                             if ring_should_be_expanded:
#                                 break
                    
#                     # 如果环需要扩展，则将所有环原子添加到待扩展集合
#                     if ring_should_be_expanded:
#                         for ring_atom_idx in ring_atom_indices:
#                             atom = mol.GetAtomWithIdx(ring_atom_idx)
#                             if atom.GetAtomMapNum():
#                                 ring_atoms_to_expand.add(atom.GetAtomMapNum())
                
#                 # 将待扩展的环原子添加到 extra_atoms 字典中
#                 for map_num in ring_atoms_to_expand:
#                     for atom in mol.GetAtoms():
#                         if atom.GetAtomMapNum() == map_num:
#                             extra_atoms[str(map_num)] = atom
#                             break
#         if expand_hetero:
#             for atom in changed_atoms:
#                 for atom2 in atom.GetNeighbors():
#                     if atom2.GetSymbol() not in ["C", "H"] and atom2.GetAtomMapNum():
#                         extra_atoms[str(atom2.GetAtomMapNum())] = atom2

#         old_tags = list(changed_atom_tags)
#         for atom_num, atom in extra_atoms.items():
#             if atom_num not in old_tags:
#                 changed_atoms.append(atom)
#                 changed_atom_tags.append(atom_num)

#         if err:
#             raise ReactionException(
#                 "Template generation failed: could not obtained changed atoms"
#             )
#         if not changed_atom_tags:
#             raise ReactionException("Template generation failed: no atoms changes")

#         # Get fragments for reactants
#         (reactant_fragments, _, _,) = extractor.get_fragments_for_changed_atoms(
#             reactants,
#             changed_atom_tags,
#             radius=radius,
#             expansion=[],
#             category="reactants",
#         )
#         # Get fragments for products %%!
#         # (WITHOUT matching groups but WITH the addition of reactant fragments)
#         product_fragments, _, _ = extractor.get_fragments_for_changed_atoms(
#             products,
#             changed_atom_tags,
#             radius=radius,
#             expansion=extractor.expand_changed_atom_tags(
#                 changed_atom_tags,
#                 reactant_fragments,
#             ),
#             category="products",
#         )

#         rxn_string = f"{reactant_fragments}>>{product_fragments}"
#         canonical_template = extractor.canonicalize_transform(rxn_string)
#         # Change from inter-molecular to intra-molecular
#         canonical_template_split = canonical_template.split(">>")
#         canonical_smarts = (
#             canonical_template_split[0][1:-1].replace(").(", ".")
#             + ">>"
#             + canonical_template_split[1][1:-1].replace(").(", ".")
#         )

#         return canonical_smarts

def smi_tokenizer(smi):
    pattern = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    assert smi == ''.join(tokens)
    tokens = [token for token in tokens if token != ':']
    return ' '.join(tokens)

def precise_tempalte_extraction(reaction_smiles, jump=3, expand_ring=True):
    template = ChemicalReaction(reaction_smiles).generate_reaction_template(jump, expand_ring=True)[0].smarts
    # template = _generate_rdchiral_template(
    #     reactants=[Chem.MolFromSmiles(smi) for smi in reaction_smiles.split(">>")[0].split(".")],
    #     products=[Chem.MolFromSmiles(smi) for smi in reaction_smiles.split(">>")[1].split(".")],
    #     radius=jump,
    #     extra_atoms={},
    #     expand_hetero=False,
    # )
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
split = "train"
jump_num = 1
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
    # conda activate unirxn
    split_num = ""
    # 50K enhance
    # base_dir = f"/root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/{split}/raw_mapping_reaction-{split}.txt"
    # cls_dir = f"/root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/{split}/reaction_cls-{split}.txt"
    # save_dir = f"/root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/{split}/tmp_smarts_{split_num}_jump{jump_num}_enhance_ring.txt"
    # save_cls_dir = f"/root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/{split}/tmp_smarts_cls_{split_num}_jump{jump_num}.txt"
    # error_dir = f"/root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/{split}/error_reaction_{split_num}_jump{jump_num}.txt"
    # 50K
    base_dir = f"/root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/{split}/raw_mapping_reaction-{split}.txt"
    cls_dir = f"/root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/{split}/reaction_cls-{split}.txt"
    save_dir = f"/root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/{split}/tmp_smarts_{split_num}_jump{jump_num}.txt"
    save_cls_dir = f"/root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/{split}/tmp_smarts_cls_{split_num}_jump{jump_num}.txt"
    error_dir = f"/root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/{split}/error_reaction_{split_num}_jump{jump_num}.txt"

    # # FULL
    # base_dir = f"/root/reaction_data/pretrain_aug/USPTO_FULL_PtoTMPtoR_aug10/{split}/raw_mapping_reaction-{split}.txt"
    # cls_dir = f"/root/reaction_data/pretrain_aug/USPTO_FULL_PtoTMPtoR_aug10/{split}/reaction_cls-{split}.txt"
    # save_dir = f"/root/reaction_data/pretrain_aug/USPTO_FULL_PtoTMPtoR_aug10/{split}/tmp_smarts_{split_num}_jump{jump_num}_enhance_ring.txt"
    # save_cls_dir = f"/root/reaction_data/pretrain_aug/USPTO_FULL_PtoTMPtoR_aug10/{split}/tmp_smarts_cls_{split_num}_jump{jump_num}.txt"
    # error_dir = f"/root/reaction_data/pretrain_aug/USPTO_FULL_PtoTMPtoR_aug10/{split}/error_reaction_{split_num}_jump{jump_num}.txt" 

    # MIT
    # base_dir = f"/root/reaction_data/pretrain_aug/USPTO_MIT_PtoTMPtoR_aug10/{split}/raw_mapping_reaction-{split}.txt"
    # cls_dir = f"/root/reaction_data/pretrain_aug/USPTO_MIT_PtoTMPtoR_aug10/{split}/reaction_cls-{split}.txt"
    # save_dir = f"/root/reaction_data/pretrain_aug/USPTO_MIT_PtoTMPtoR_aug10/{split}/tmp_smarts_{split_num}_jump{jump_num}_enhance_ring.txt"
    # save_cls_dir = f"/root/reaction_data/pretrain_aug/USPTO_MIT_PtoTMPtoR_aug10/{split}/tmp_smarts_cls_{split_num}_jump{jump_num}.txt"
    # error_dir = f"/root/reaction_data/pretrain_aug/USPTO_MIT_PtoTMPtoR_aug10/{split}/error_reaction_{split_num}_jump{jump_num}.txt"  

    # Read input data
    with open(base_dir, 'r') as f:
        data = [line.strip() for line in f.readlines() if line.strip()]
    # with open(cls_dir, 'r') as f:
    #     cls_data = [line.strip() for line in f.readlines() if line.strip()]
    cls_data = ["999"] * len(data)
    assert len(data) == len(cls_data), "The number of reactions must match the number of classes."
    # Initialize multiprocessing pool with 40 cores
    pool = multiprocessing.Pool(processes=45)
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
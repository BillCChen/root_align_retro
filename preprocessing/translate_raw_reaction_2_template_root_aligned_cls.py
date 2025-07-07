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


split = "train"
base_dir = f"/root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/{split}/tmp_smarts.txt"
cls_dir = f"/root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/{split}/reaction_cls-{split}.txt"
# 添加的代码：创建分类文件夹
cls_base_dir = f"/root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/{split}/cls"
os.makedirs(cls_base_dir, exist_ok=True)
for cls_id in range(1, 11):
    os.makedirs(f"{cls_base_dir}/cls_{cls_id}", exist_ok=True)

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
import sys
sys.path.append("/root/retro_synthesis/template_analysis")
from tools.handle_templates import Get_Reaction_
# 添加 /root/retro_synthesis/reaction_utils 的库查找路径
sys.path.append("/root/retro_synthesis/reaction_utils")
from rxnutils.chem.reaction import ChemicalReaction
def smi_tokenizer(smi):
    pattern = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    assert smi == ''.join(tokens)
    return ' '.join(tokens)

def precise_tempalte_extraction(reaction_smiles,jump=3,expand_ring=True):
    template =  ChemicalReaction(reaction_smiles).generate_reaction_template(jump,expand_ring)[0].smarts
    # 确保不要浪费内存
    assert type(template) == str, "Template extraction failed, expected a string."
    return template


def clear_map_canonical_smiles(smi, canonical=True, root=-1, type_='smarts'):
    assert type_ in ['smiles', 'smarts'], "Type must be either 'smiles' or 'smarts'."
    if type_ == 'smiles':
        mol = Chem.MolFromSmiles(smi,sanitize=False)
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
    """
    Remove atom mapping from a reaction template.
    Args:
        template (str): Reaction template with atom mapping.
        mode: 
            single: only use function_clear_map_canonical_smiles
            double: use function_clear_map_canonical_smiles and Chem.MolToSmiles(Chem.MolFromSmiles()
    Returns:
        str: Reaction template without atom mapping.
    """
    def remove_mapping(smi):
        return clear_map_canonical_smiles(smi, canonical=True, type_='smarts')

    left , right = template.split('>>')
    if len(left.split('.')) > 1:
        left = '.'.join([remove_mapping(c) for c in left.split('.')])
    else:
        left = remove_mapping(left)
    if len(right.split('.')) > 1:
        right = '.'.join([remove_mapping(c) for c in right.split('.')])
    else:
        right = remove_mapping(right)
    return left + '>>' + right
def get_root_id(mol, root_map_number):
    """获取指定原子映射号的原子索引"""
    for i, atom in enumerate(mol.GetAtoms()):
        if atom.GetAtomMapNum() == root_map_number:
            return i
    return -1

def get_cano_map_number(smi, root=-1,type_='smarts'):
    """获取原子映射号的规范顺序（用于对齐）"""
    assert type_ in ['smiles', 'smarts'], "Type must be either 'smiles' or 'smarts'."
    if type_ == 'smiles':
        atommap_mol = Chem.MolFromSmiles(smi)
        canonical_mol = Chem.MolFromSmiles(clear_map_canonical_smiles(smi, root=root))
    elif type_ == 'smarts':
        atommap_mol = Chem.MolFromSmarts(smi)
        canonical_mol = Chem.MolFromSmarts(clear_map_canonical_smiles(smi, root=root, type_='smarts'))
    cano2atommapIdx = atommap_mol.GetSubstructMatch(canonical_mol)
    if len(cano2atommapIdx) < canonical_mol.GetNumAtoms():
        return None
    return [atommap_mol.GetAtomWithIdx(idx).GetAtomMapNum() for idx in cano2atommapIdx]

# 修改后的模板生成逻辑
def generate_root_aligned_templates(reactant_smi, product_smi, augmentation=1):
    """生成 root-aligned 的反应模板"""
    # pro_mol = Chem.MolFromSmiles(product_smi, sanitize=False)
    # rea_mols = [Chem.MolFromSmiles(smi, sanitize=False) for smi in reactant_smi.split('.')]
    pro_mol = Chem.MolFromSmarts(product_smi)
    rea_mols = [Chem.MolFromSmarts(smi) for smi in reactant_smi.split('.')]
    
    # 获取产物原子的映射号（用于选择不同的根原子）
    pro_atom_map_numbers = list(map(int, re.findall(r"(?<=:)\d+", product_smi)))
    
    templates = []
    for k in range(augmentation):
        try:
            # 随机选择一个产物原子作为根（或按顺序选择）
            pro_root_map = random.choice(pro_atom_map_numbers) if augmentation > 1 else pro_atom_map_numbers[0]
            pro_root = get_root_id(pro_mol, pro_root_map)
            
            # 生成 root-aligned 的产物 SMILES（无原子映射）
            aligned_product = clear_map_canonical_smiles(product_smi, root=pro_root, type_='smarts')
            
            # 对齐反应物（基于产物原子的映射）
            aligned_reactants = []
            cano_atom_map = get_cano_map_number(product_smi, root=pro_root)
            if cano_atom_map is None:
                continue  # 跳过无法对齐的情况
            
            for rea_smi in reactant_smi.split('.'):
                rea_map_numbers = list(map(int, re.findall(r"(?<=:)\d+", rea_smi)))
                for map_num in cano_atom_map:
                    if map_num in rea_map_numbers:
                        # rea_root = get_root_id(Chem.MolFromSmiles(rea_smi), map_num)
                        rea_root = get_root_id(Chem.MolFromSmarts(rea_smi), map_num)
                        aligned_rea = clear_map_canonical_smiles(rea_smi, root=rea_root, type_='smarts')
                        aligned_reactants.append(aligned_rea)
                        break

            # 组合成反应模板
            if aligned_reactants:
                template = f"{'.'.join(aligned_reactants)}>>{aligned_product}"
                templates.append(template)

        except Exception as e:
            print(f"Using original {k} reaction {reactant_smi[:10]} >> {product_smi[:10]} as fallback.")
            # templates.append(f"{reactant_smi}>>{product_smi}")  # 保留原始反应作为备选
    
    return templates

# 主流程修改
with open(base_dir, 'r') as f:
    data = f.readlines()
data = [line.strip() for line in data if line.strip()]

# 读取分类信息
with open(cls_dir, 'r') as f:
    cls_data = f.readlines()
cls_data = [line.strip() for line in cls_data if line.strip()]
# 每 20个数据取一个
cls_data = cls_data[::20]
assert len(cls_data) == len(data), "The number of classification data must match the number of reactions."

from collections import defaultdict
# 使用字典按类别存储模板
cls_templates = defaultdict(lambda: {'product': [], 'reactant': []})

for reaction, cls_label in tqdm(zip(data, cls_data), desc="Processing reactions by class", total=len(data)):
    try:
        cls_id = int(cls_label)
        reactant, product = reaction.split('>>')
        
        templates = generate_root_aligned_templates(reactant, product, augmentation=20)
        if len(templates) == 0:
            print(f"Warning: No templates generated for reaction: {reaction[:10]}, Class: {cls_label}")
            continue
        if len(templates) != 20:
            print(f"Warning: Expected 20 templates, but got {len(templates)} for reaction: {reaction}")
            # 随机复制已有的进行补全
        while len(templates) < 20:
            templates.append(random.choice(templates))
        for template in templates:
            tmp_reactant, tmp_product = template.split('>>')
            # Token化
            token_reactant = smi_tokenizer(tmp_reactant)
            token_product = smi_tokenizer(tmp_product)
            
            cls_templates[cls_id]['reactant'].append(token_reactant)
            cls_templates[cls_id]['product'].append(token_product)
            
    except Exception as e:
        print(f"Error processing reaction: {reaction}, Class: {cls_label}, Error: {e}")

# 保存分类数据
for cls_id, templates in cls_templates.items():
    cls_dir = f"{cls_base_dir}/cls_{cls_id}"
    
    with open(f"{cls_dir}/tmp_product.txt", "w") as f:
        f.write("\n".join(templates['product']))
    
    with open(f"{cls_dir}/tmp_reactant.txt", "w") as f:
        f.write("\n".join(templates['reactant']))

print("数据处理完成，已按类别保存")
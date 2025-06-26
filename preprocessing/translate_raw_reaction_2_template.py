split = "train"
base_dir = f"/root/retro_synthesis/root_align_retro/dataset/USPTO_50K_PtoTMPtoR_aug20/{split}/raw_mapping_reaction-{split}.txt"
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
with open(base_dir, 'r') as f:
    data = f.readlines()
data = [line.strip() for line in data if line.strip()]
all_tmp_product = []
all_tmp_reactant = []
store_temp = ""
from tqdm import tqdm
for reaction_ in tqdm(data, desc="Processing reactions", ncols=100):
    reactant , product, num = reaction_.split('>>')
    if num == "0":
        # print(f"{reactant}>>{product}")
        store_temp = remove_atom_mapping_template(precise_tempalte_extraction(f"{reactant}>>{product}",3))
        tmp_reactant,tmp_prodcut = store_temp.split('>>')
        assert len(tmp_reactant) > 0, "The reactant part of the template cannot be empty."
        assert len(tmp_prodcut) > 0, "The product part of the template cannot be empty."
        all_tmp_reactant.append(tmp_reactant)
        all_tmp_product.append(tmp_prodcut)
        # print(f"{tmp_reactant}>>{tmp_prodcut}")
    else:
        tmp_reactant, tmp_prodcut = store_temp.split('>>')
        all_tmp_reactant.append(tmp_reactant)
        all_tmp_product.append(tmp_prodcut)
assert len(all_tmp_reactant) == len(all_tmp_product), "The number of reactants and products must be the same."
assert len(all_tmp_reactant) == len(data), "The number of reactants and products must match the number of reactions."
assert len(all_tmp_product) == len(data), "The number of reactants and products must match the number of reactions."
all_tmp_product = [smi_tokenizer(smi) for smi in all_tmp_product]
all_tmp_reactant = [smi_tokenizer(smi) for smi in all_tmp_reactant]
with open(f"/root/retro_synthesis/root_align_retro/dataset/USPTO_50K_PtoTMPtoR_aug20/{split}/tmp_product.txt", "w") as f:
    for smi in all_tmp_product:
        f.write(smi + '\n')
with open(f"/root/retro_synthesis/root_align_retro/dataset/USPTO_50K_PtoTMPtoR_aug20/{split}/tmp_reactant.txt", "w") as f:
    for smi in all_tmp_reactant:
        f.write(smi + '\n')
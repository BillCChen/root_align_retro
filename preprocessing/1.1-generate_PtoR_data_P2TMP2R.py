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
# import sys
# sys.path.append("/root/retro_synthesis/template_analysis")
# from tools.handle_templates import Get_Reaction_
# # 添加 /root/retro_synthesis/reaction_utils 的库查找路径
# sys.path.append("/root/retro_synthesis/reaction_utils")
# from rxnutils.chem.reaction import ChemicalReaction
def smi_tokenizer(smi):
    pattern = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    assert smi == ''.join(tokens)
    return ' '.join(tokens)

# def precise_tempalte_extraction(reaction_smiles,jump=3,expand_ring=True):
#     template =  ChemicalReaction(reaction_smiles).generate_reaction_template(jump,expand_ring)[0].smarts
#     # 确保不要浪费内存
#     assert type(template) == str, "Template extraction failed, expected a string."
#     return template


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
def remove_atom_mapping_template(template,mode="double"):
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
        if mode == "single":
            return clear_map_canonical_smiles(smi, canonical=True, type_='smiles')
        elif mode == "double":
            return Chem.MolToSmiles(Chem.MolFromSmiles(clear_map_canonical_smiles(smi, canonical=True, type_='smiles'), sanitize=False),rootedAtAtom=-1, canonical=True)
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

def get_cano_map_number(smi,root=-1):
    atommap_mol = Chem.MolFromSmiles(smi)
    canonical_mol = Chem.MolFromSmiles(clear_map_canonical_smiles(smi,root=root))
    cano2atommapIdx = atommap_mol.GetSubstructMatch(canonical_mol)
    correct_mapped = [canonical_mol.GetAtomWithIdx(i).GetSymbol() == atommap_mol.GetAtomWithIdx(index).GetSymbol() for i,index in enumerate(cano2atommapIdx)]
    atom_number = len(canonical_mol.GetAtoms())
    if np.sum(correct_mapped) < atom_number or len(cano2atommapIdx) < atom_number:
        cano2atommapIdx = [0] * atom_number
        atommap2canoIdx = canonical_mol.GetSubstructMatch(atommap_mol)
        if len(atommap2canoIdx) != atom_number:
            return None
        for i, index in enumerate(atommap2canoIdx):
            cano2atommapIdx[index] = i
    id2atommap = [atom.GetAtomMapNum() for atom in atommap_mol.GetAtoms()]

    return [id2atommap[cano2atommapIdx[i]] for i in range(atom_number)]

def get_root_id(mol,root_map_number):
    root = -1
    for i, atom in enumerate(mol.GetAtoms()):
        if atom.GetAtomMapNum() == root_map_number:
            root = i
            break
    return root
    # root = -1
    # for i, atom in enumerate(mol.GetAtoms()):
    #     if atom.GetAtomMapNum() == root_map_number:
    #         return i

"""multiprocess"""
def preprocess(save_dir, reactants, products,set_name, augmentation=1, reaction_types=None,root_aligned=True,character=False, processes=-1):
    """
    preprocess reaction data to extract graph adjacency matrix and features
    """

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    data = [ {
        "reactant":i,
        "product":j,
        "augmentation":augmentation,
        "root_aligned":root_aligned,
    }  for i,j in zip(reactants,products)]
    src_data = []
    tgt_data = []
    # tmp_data = []
    # tmp_data2 = []
    raw_mapping_data = []
    skip_dict = {
        'invalid_p':0,
        'invalid_r':0,
        'small_p':0,
        'small_r':0,
        'error_mapping':0,
        'error_mapping_p':0,
        'empty_p':0,
        'empty_r':0,
        'error_main_function':0
    }
    processes = multiprocessing.cpu_count() if processes < 0 else processes
    pool = multiprocessing.Pool(processes=processes)
    results = pool.map(func=multi_process,iterable=data)
    pool.close()
    pool.join()
    edit_distances = []
    edit_distances_std = []
    for result in tqdm(results):
        if result['status'] != 0:
            skip_dict[result['status']] += 1
            continue
        if character:
            for i in range(len(result['src_data'])):
                result['src_data'][i] = " ".join([char for char in "".join(result['src_data'][i].split())])
            for i in range(len(result['tgt_data'])):
                result['tgt_data'][i] = " ".join([char for char in "".join(result['tgt_data'][i].split())])
            # for i in range(len(result['tmp_data'])):
            #     result['tmp_data'][i] = " ".join([char for char in "".join(result['tmp_data'][i].split())])
            # for i in range(len(result['tmp_data2'])):
            #     result['tmp_data2'][i] = " ".join([char for char in "".join(result['tmp_data2'][i].split())])
        edit_distances.append(result['edit_distance'])
        edit_distances_std.append(result['edit_distance_std'])
        src_data.extend(result['src_data'])
        tgt_data.extend(result['tgt_data'])
        # tmp_data.extend(result['tmp_data'])
        # tmp_data2.extend(result['tmp_data2'])
        raw_mapping_data.extend(result['raw_mapping_data'])

    print("Avg. edit distance:", np.mean(edit_distances))
    print("Std. edit distance:", np.mean(edit_distances_std))
    print('size', len(src_data))
    for key,value in skip_dict.items():
        print(f"{key}:{value},{value/len(reactants)}")
    if augmentation != 999:
        with open(
                os.path.join(save_dir, 'src-{}.txt'.format(set_name)), 'w') as f:
            for src in src_data:
                f.write('{}\n'.format(src))

        with open(
                os.path.join(save_dir, 'tgt-{}.txt'.format(set_name)), 'w') as f:
            for tgt in tgt_data:
                f.write('{}\n'.format(tgt))
        # with open(
        #         os.path.join(save_dir, 'tmp-{}.txt'.format(set_name)), 'w') as f:
        #     for tmp in tmp_data:
        #         f.write('{}\n'.format(tmp))
        # with open(
        #         os.path.join(save_dir, 'tmp2-{}.txt'.format(set_name)), 'w') as f:
        #     for tmp2 in tmp_data2:
        #         f.write('{}\n'.format(tmp2))
        with open(
                os.path.join(save_dir, 'raw_mapping_reaction-{}.txt'.format(set_name)), 'w') as f:
            for raw_mapping in raw_mapping_data:
                f.write('{}\n'.format(raw_mapping))
    return src_data,tgt_data # tmp_data, tmp_data2


def multi_process(data):

    pt = re.compile(r':(\d+)]')
    product = data['product']
    reactant = data['reactant']
    augmentation = data['augmentation']
    # pro_mol = Chem.MolFromSmiles(product)
    # rea_mol = Chem.MolFromSmiles(reactant)
    pro_mol = Chem.MolFromSmarts(product)
    rea_mol = Chem.MolFromSmarts(reactant)
    """checking data quality"""
    rids = sorted(re.findall(pt, reactant))
    pids = sorted(re.findall(pt, product))
    return_status = {
        "status":0,
        "src_data":[],
        "tgt_data":[],
        # "tmp_data":[],
        # "tmp_data2":[],
        'raw_mapping_data':[],
        "edit_distance":0,
        "edit_distance_std":0,
    }
    if ",".join(rids) != ",".join(pids):  # mapping is not 1:1
        return_status["status"] = "error_mapping"
    if len(set(rids)) != len(rids):  # mapping is not 1:1
        return_status["status"] = "error_mapping"
    if len(set(pids)) != len(pids):  # mapping is not 1:1
        return_status["status"] = "error_mapping"
    if "" == product:
        return_status["status"] = "empty_p"
    try:
        if len(rea_mol.GetAtoms()) < 5:
            return_status["status"] = "invalid_r"
    except:
        return_status["status"] = "invalid_r"
    try:
        if len(pro_mol.GetAtoms()) < 5:
            return_status["status"] = "invalid_p"
    except:
        return_status["status"] = "invalid_p"
    if "" == reactant:
        return_status["status"] = "empty_r"
    if rea_mol is None:
        return_status["status"] = "invalid_r"
    # if len(rea_mol.GetAtoms()) < 5:
        # return_status["status"] = "small_r"
    if pro_mol is None:
        return_status["status"] = "invalid_p"
    # if len(pro_mol.GetAtoms()) == 1:
    #     return_status["status"] = "small_p"
    if not all([a.HasProp('molAtomMapNumber') for a in pro_mol.GetAtoms()]):
        return_status["status"] = "error_mapping_p"
    """finishing checking data quality"""
    try:
        if return_status['status'] == 0:
            pro_atom_map_numbers = list(map(int, re.findall(r"(?<=:)\d+", product)))
            reactant = reactant.split(".")
            if data['root_aligned']:
                reversable = False  # no shuffle
                # augmentation = 100
                if augmentation == 999:
                    product_roots = pro_atom_map_numbers
                    times = len(product_roots)
                else:
                    product_roots = [-1]
                    # reversable = len(reactant) > 1

                    max_times = len(pro_atom_map_numbers)
                    times = min(augmentation, max_times)
                    if times < augmentation:  # times = max_times
                        product_roots.extend(pro_atom_map_numbers)
                        product_roots.extend(random.choices(product_roots, k=augmentation - len(product_roots)))
                    else:  # times = augmentation
                        while len(product_roots) < times:
                            product_roots.append(random.sample(pro_atom_map_numbers, 1)[0])
                            # pro_atom_map_numbers.remove(product_roots[-1])
                            if product_roots[-1] in product_roots[:-1]:
                                product_roots.pop()
                    times = len(product_roots)
                    assert times == augmentation
                    if reversable:
                        times = int(times / 2)
                # candidates = []
                # template = remove_atom_mapping_template(precise_tempalte_extraction(f"{'.'.join(reactant)}>>{product}",3))
                
                # tmp2,tmp = template.split('>>')
                for k in range(times):
                    pro_root_atom_map = product_roots[k]
                    pro_root = get_root_id(pro_mol, root_map_number=pro_root_atom_map)
                    cano_atom_map = get_cano_map_number(product, root=pro_root)
                    if cano_atom_map is None:
                        return_status["status"] = "error_mapping"
                        return return_status
                    pro_smi = clear_map_canonical_smiles(product, canonical=True, root=pro_root)
                    aligned_reactants = []
                    aligned_reactants_order = []
                    rea_atom_map_numbers = [list(map(int, re.findall(r"(?<=:)\d+", rea))) for rea in reactant]
                    used_indices = []
                    for i, rea_map_number in enumerate(rea_atom_map_numbers):
                        for j, map_number in enumerate(cano_atom_map):
                            # select mapping reactans
                            if map_number in rea_map_number:
                                # rea_root = get_root_id(Chem.MolFromSmiles(reactant[i]), root_map_number=map_number)
                                rea_root = get_root_id(Chem.MolFromSmarts(reactant[i]), root_map_number=map_number)
                                rea_smi = clear_map_canonical_smiles(reactant[i], canonical=True, root=rea_root)
                                aligned_reactants.append(rea_smi)
                                aligned_reactants_order.append(j)
                                used_indices.append(i)
                                break
                    sorted_reactants = sorted(list(zip(aligned_reactants, aligned_reactants_order)), key=lambda x: x[1])
                    aligned_reactants = [item[0] for item in sorted_reactants]
                    reactant_smi = ".".join(aligned_reactants)
                    product_tokens = smi_tokenizer(pro_smi)
                    # print("product:", pro_smi)
                    # print("product:", product_tokens)
                    reactant_tokens = smi_tokenizer(reactant_smi)
                    # tmp_tokens = smi_tokenizer(tmp)
                    # tmp2_tokens = smi_tokenizer(tmp2)
                    # if len(tmp_tokens) > 1 and len(tmp2_tokens) > 1:
                    return_status['src_data'].append(product_tokens)
                    return_status['tgt_data'].append(reactant_tokens)
                    return_status['raw_mapping_data'].append(f"{'.'.join(reactant)}>>{product}>>{k}")
                        # return_status['tmp_data'].append(tmp_tokens)
                        # return_status['tmp_data2'].append(tmp2_tokens)
                    # else:
                    #     print(f"Original: {pro_smi}>>{reactant_smi}")
                    #     print(f"Template {tmp} or {tmp2} is too short, skipping.")
                    if reversable:
                        aligned_reactants.reverse()
                        reactant_smi = ".".join(aligned_reactants)
                        product_tokens = smi_tokenizer(pro_smi)
                        reactant_tokens = smi_tokenizer(reactant_smi)
                        return_status['src_data'].append(product_tokens)
                        return_status['tgt_data'].append(reactant_tokens)
                        return_status['raw_mapping_data'].append(f"{'.'.join(reactant)}>>{product}>>{k}")
                assert len(return_status['src_data']) == data['augmentation']
            else:
                cano_product = clear_map_canonical_smiles(product)
                cano_reactanct = ".".join([clear_map_canonical_smiles(rea) for rea in reactant if len(set(map(int, re.findall(r"(?<=:)\d+", rea))) & set(pro_atom_map_numbers)) > 0 ])
                return_status['src_data'].append(smi_tokenizer(cano_product))
                return_status['tgt_data'].append(smi_tokenizer(cano_reactanct))
                pro_mol = Chem.MolFromSmiles(cano_product)
                rea_mols = [Chem.MolFromSmiles(rea) for rea in cano_reactanct.split(".")]
                for i in range(int(augmentation-1)):
                    pro_smi = Chem.MolToSmiles(pro_mol,doRandom=True)
                    rea_smi = [Chem.MolToSmiles(rea_mol,doRandom=True) for rea_mol in rea_mols]
                    rea_smi = ".".join(rea_smi)
                    return_status['src_data'].append(smi_tokenizer(pro_smi))
                    return_status['tgt_data'].append(smi_tokenizer(rea_smi))
            edit_distances = []
            for src,tgt in zip(return_status['src_data'],return_status['tgt_data']):
                edit_distances.append(textdistance.levenshtein.distance(src.split(),tgt.split()))
            return_status['edit_distance'] = np.mean(edit_distances)
            return_status['edit_distance_std'] = np.std(edit_distances)
        
    except:
        return_status["status"] = "error_main_function"
        
    
    return return_status


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-dataset',
                        type=str,
                        default='USPTO_50K')
    parser.add_argument("-augmentation",type=int,default=20)
    parser.add_argument("-seed",type=int,default=33)
    parser.add_argument("-processes",type=int,default=20)
    parser.add_argument("-test_only", action="store_true")
    parser.add_argument("-train_only", action="store_true")
    parser.add_argument("-test_except", action="store_true")
    parser.add_argument("-validastrain", action="store_true")
    parser.add_argument("-character", action="store_true")
    parser.add_argument("-canonical", action="store_true")
    parser.add_argument("-postfix",type=str,default="")
    args = parser.parse_args()
    print('preprocessing dataset {}...'.format(args.dataset))
    # assert args.dataset in ['USPTO_50K', 'USPTO_full','USPTO-MIT']
    print(args)
    if args.test_only:
        datasets = ['test']
    elif args.train_only:
        datasets = ['train']
    elif args.test_except:
        datasets = ['val', 'train']
    elif args.validastrain:
        datasets = ['test', 'val', 'train']
    else:
        datasets = ['test', 'val', 'train']

    random.seed(args.seed)
    if args.dataset == "USPTO-FULL":
        datadir = '/root/reaction_data/pretrain_aug/USPTO_FULL_ori/{}'.format(args.dataset)
        savedir = '/root/reaction_data/pretrain_aug/USPTO_FULL_ori/{}_PtoR_aug{}'.format(args.dataset,args.augmentation)
        # datadir = "/root/reaction_data/pretrain_aug/USPTO_full_ori_txt"
        # savedir = '/root/reaction_data/pretrain_aug/USPTO_full_PtoTMPtoR_aug{}'.format(args.augmentation)
        savedir += args.postfix
        if not os.path.exists(savedir):
            os.makedirs(savedir)
        for i, data_set in enumerate(datasets):
            with open(os.path.join(datadir,f"{data_set}.txt"),"r") as f:
                reaction_list = f.readlines()
                if args.validastrain and data_set == "train":
                    with open(os.path.join(datadir, f"valid.txt"), "r") as f:
                        reaction_list += f.readlines()

                reactant_smarts_list = list(
                    map(lambda x: x.split('>>')[0], reaction_list))
                product_smarts_list = list(
                    map(lambda x: x.split('>>')[1], reaction_list))
                product_smarts_list = list(
                    map(lambda x: x.split(' ')[0], product_smarts_list))
                save_dir = os.path.join(savedir, data_set)

                multiple_product_indices = [i for i in range(len(product_smarts_list)) if "." in product_smarts_list[i]]
                for index in multiple_product_indices:
                    products = product_smarts_list[index].split(".")
                    for product in products:
                        reactant_smarts_list.append(reactant_smarts_list[index])
                        product_smarts_list.append(product)
                for index in multiple_product_indices[::-1]:
                    del reactant_smarts_list[index]
                    del product_smarts_list[index]

                src_data, tgt_data = preprocess(
                    save_dir,
                    reactant_smarts_list,
                    product_smarts_list,
                    data_set,
                    args.augmentation,
                    reaction_types=None,
                    root_aligned=not args.canonical,
                    character=args.character,
                    processes=args.processes,
                )
    elif args.dataset == "CHORISO":
        datadir = '/root/reaction_data/pretrain_aug/CHORISO_PtoTMPtoR_aug10/ori_data'
        savedir = '/root/reaction_data/pretrain_aug/CHORISO_PtoTMPtoR_aug{}/'.format(args.augmentation)
        savedir += args.postfix
        if not os.path.exists(savedir):
            os.makedirs(savedir)
        for i, data_set in enumerate(datasets):
            with open(os.path.join(datadir,f"{data_set}.txt"),"r") as f:
                reaction_list = f.readlines()
                if args.validastrain and data_set == "train":
                    with open(os.path.join(datadir, f"valid.txt"), "r") as f:
                        reaction_list += f.readlines()
                # reaction_list = reaction_list[:10000]
                reactant_smarts_list = list(
                    map(lambda x: x.split('>>')[0], reaction_list))
                product_smarts_list = list(
                    map(lambda x: x.split('>>')[1], reaction_list))
                product_smarts_list = list(
                    map(lambda x: x.split(' ')[0], product_smarts_list))
                save_dir = os.path.join(savedir, data_set)
                print("Total Data Size", len(reaction_list))
                multiple_product_indices = [i for i in range(len(product_smarts_list)) if "." in product_smarts_list[i]]
                for index in multiple_product_indices:
                    products = product_smarts_list[index].split(".")
                    for product in products:
                        reactant_smarts_list.append(reactant_smarts_list[index])
                        product_smarts_list.append(product)
                for index in multiple_product_indices[::-1]:
                    del reactant_smarts_list[index]
                    del product_smarts_list[index]

                src_data, tgt_data = preprocess(
                    save_dir,
                    reactant_smarts_list,
                    product_smarts_list,
                    data_set,
                    args.augmentation,
                    reaction_types=None,
                    root_aligned=not args.canonical,
                    character=args.character,
                    processes=args.processes,
                )
    else:
        datadir = '/root/reaction_data/pretrain_aug/{}'.format(args.dataset)
        savedir = '/root/reaction_data/pretrain_aug/{}_PtoTMPtoR_aug{}_trash'.format(args.dataset, args.augmentation)

        savedir += args.postfix
        if not os.path.exists(savedir):
            os.makedirs(savedir)
        for i, data_set in enumerate(datasets):
            csv_path = f"{datadir}/raw_{data_set}.csv"
            csv = pd.read_csv(csv_path)
            reaction_list = list(csv["reactants>reagents>production"])
            if args.validastrain and data_set == "train":
                csv_path = f"{datadir}/raw_val.csv"
                csv = pd.read_csv(csv_path)
                reaction_list += list(csv["reactants>reagents>production"])

            # random.shuffle(reaction_list)
            reactant_smarts_list = list(
                map(lambda x: x.split('>')[0], reaction_list))
            reactant_smarts_list = list(
                map(lambda x: x.split(' ')[0], reactant_smarts_list))
            reagent_smarts_list = list(
                map(lambda x: x.split('>')[1], reaction_list))
            product_smarts_list = list(
                map(lambda x: x.split('>')[2], reaction_list))
            product_smarts_list = list(
                map(lambda x: x.split(' ')[0], product_smarts_list))  # remove ' |f:1...'
            print("Total Data Size", len(reaction_list))

            # reaction_class_list = list(map(lambda x: int(x) - 1, csv['class']))
            sub_react_list = reactant_smarts_list
            sub_prod_list = product_smarts_list
            save_dir = os.path.join(savedir, data_set)

            # duplicate multiple product reactions into multiple ones with one product each
            multiple_product_indices = [i for i in range(len(sub_prod_list)) if "." in sub_prod_list[i]]
            for index in multiple_product_indices:
                products = sub_prod_list[index].split(".")
                for product in products:
                    sub_react_list.append(sub_react_list[index])
                    sub_prod_list.append(product)
            for index in multiple_product_indices[::-1]:
                del sub_react_list[index]
                del sub_prod_list[index]
            src_data, tgt_data = preprocess(
                save_dir,
                sub_react_list,
                sub_prod_list,
                data_set,
                args.augmentation,
                reaction_types=None,
                root_aligned=not args.canonical,
                character=args.character,
                processes=args.processes,
            )

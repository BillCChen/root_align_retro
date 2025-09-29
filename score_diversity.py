
import os
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "1"
import re
import torch
from onmt.translate.translator import build_translator
from types import SimpleNamespace
import argparse
import os
import logging
from rdkit import Chem
class RetroSynthesisDemo:
    def __init__(self, model_path, gpu_device=0, beam_size=10, n_best=3, batch_size=24):
        self.model_path = model_path
        self.gpu_device = gpu_device
        self.beam_size = beam_size
        self.n_best = n_best
        self.batch_size = batch_size

        # 构建 OpenNMT 参数
        opt = SimpleNamespace(
            models=[self.model_path],
            gpu=self.gpu_device,
            beam_size=self.beam_size,
            n_best=self.n_best,
            batch_size=self.batch_size,
            max_length=256,
            seed=42,
            # === 需要补充的字段 ===
            block_ngram_repeat=0,
            ignore_when_blocking=[],
            replace_unk=True,
            verbose=False,
            report_align=False,
            report_time=False,
            attn_debug=False,
            align_debug=False,
            dump_beam="",
            ban_unk_token=False,
            phrase_table="",
            log_file="",
            log_file_level="0",
            batch_type="sents",
            min_length=0,
            max_sent_length=None,
            coverage_penalty="none",
            alpha=0.0,
            beta=0.0,
            stepwise_penalty=False,
            length_penalty="none",
            ratio=0.0,
            random_sampling_topk=0,
            random_sampling_topp=0.0,
            random_sampling_temp=1.0,
            avg_raw_probs=False,
            data_type="text",
            src=None,
            src_feats=None,
            tgt=None,
            tgt_prefix=False,
            shard_size=0,
            output="/root/z-trash/onmt_out.txt",
            fp32=False,
            int8=False,
        )

        # 构建 Translator（只加载一次）
    

        # 构建 Translator（只加载一次）
        self.inference_model = build_translator(opt, report_score=False, out_file=open(os.devnull, "w"))

    def smi_tokenizer(self,smi):
        pattern = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
        regex = re.compile(pattern)
        tokens = [token for token in regex.findall(smi)]
        # 删去 : 这个token
        tokens = [token for token in tokens if token != ':']
        assert smi == ''.join(tokens)
        return ' '.join(tokens)
    def inference(self, input_smiles):
        # 对输入做 token 化
        tokenized = [self.smi_tokenizer(smi) for smi in input_smiles]

        # 模型推理
        scores, predictions = self.inference_model.translate(
            src=tokenized,
            tgt=None,
            batch_size=self.batch_size,
            attn_debug=False
        )

        # 整理输出
        translated = []
        for pred_list in predictions:  # 每个输入样本
            for smi in pred_list:      # n_best 个预测
                translated.append(smi.replace(" ", ""))

        return translated
class ForwardSynthesisDemo:
    def __init__(self, model_path, gpu_device=0, beam_size=10, n_best=1, batch_size=24):
        self.model_path = model_path
        self.gpu_device = gpu_device
        self.beam_size = beam_size
        self.n_best = n_best
        self.batch_size = batch_size

        # 构建 OpenNMT 参数
        opt = SimpleNamespace(
            models=[self.model_path],
            gpu=self.gpu_device,
            beam_size=self.beam_size,
            n_best=self.n_best,
            batch_size=self.batch_size,
            max_length=256,
            seed=42,
            # === 需要补充的字段 ===
            block_ngram_repeat=0,
            ignore_when_blocking=[],
            replace_unk=True,
            verbose=False,
            report_align=False,
            report_time=False,
            attn_debug=False,
            align_debug=False,
            dump_beam="",
            ban_unk_token=False,
            phrase_table="",
            log_file="",
            log_file_level="0",
            batch_type="sents",
            min_length=0,
            max_sent_length=None,
            coverage_penalty="none",
            alpha=0.0,
            beta=0.0,
            stepwise_penalty=False,
            length_penalty="none",
            ratio=0.0,
            random_sampling_topk=0,
            random_sampling_topp=0.0,
            random_sampling_temp=1.0,
            avg_raw_probs=False,
            data_type="text",
            src=None,
            src_feats=None,
            tgt=None,
            tgt_prefix=False,
            shard_size=0,
            output="/root/z-trash/onmt_out.txt",
            fp32=False,
            int8=False,
        )

        # 构建 Translator（只加载一次）
    

        # 构建 Translator（只加载一次）
        self.inference_model = build_translator(opt, report_score=False, out_file=open(os.devnull, "w"))

    def smi_tokenizer(self,smi):
        pattern = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
        regex = re.compile(pattern)
        tokens = [token for token in regex.findall(smi)]
        # 删去 : 这个token
        tokens = [token for token in tokens if token != ':']
        assert smi == ''.join(tokens)
        return ' '.join(tokens)
    def inference(self, input_smiles):
        # 对输入做 token 化
        tokenized = [self.smi_tokenizer(smi) for smi in input_smiles]

        # 模型推理
        scores, predictions = self.inference_model.translate(
            src=tokenized,
            tgt=None,
            batch_size=self.batch_size,
            attn_debug=False
        )

        # 整理输出
        translated = []
        for pred_list in predictions:  # 每个输入样本
            for smi in pred_list:      # n_best 个预测
                translated.append(smi.replace(" ", ""))

        return translated

from tqdm import tqdm, trange
import sys
from easydict import EasyDict as edict
import yaml
from yaml import Dumper, Loader
cfg = edict({
    'predictor':
    yaml.load(open('/root/retro_synthesis/Uni-RXN-official-backup/config/predictor/model.yaml'),
                Loader=Loader),
})
sys.path.append('/root/retro_synthesis/Uni-RXN-official-backup/LocalTransform')
from Synthesis import init_LocalTransform, predict_product_batch
import random
from collections import defaultdict, deque
def random_substructure(smiles, r=4, d=2):
    """
    输入一个RDKit分子对象，拓展半径r和搜索半径d，随机输出一个子结构的SMILES字符串。
    
    参数:
    mol: RDKit的Mol对象
    r: 整数，拓展半径
    d: 整数，搜索半径
    
    返回:
    子结构的SMILES字符串
    """
    mol = Chem.MolFromSmiles(smiles)
    # 如果分子中没有键，则返回整个分子的SMILES
    if mol.GetNumBonds() == 0:
        return Chem.MolToSmiles(mol)
    
    # 预计算所有芳香环信息
    rings = mol.GetRingInfo().AtomRings()
    aromatic_rings = []
    for ring in rings:
        is_aromatic = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if not atom.GetIsAromatic():
                is_aromatic = False
                break
        if is_aromatic:
            aromatic_rings.append(ring)
    
    # 构建原子到所属芳香环所有原子的映射
    atom_to_aromatic_atoms = defaultdict(set)
    for ring in aromatic_rings:
        for idx1 in ring:
            for idx2 in ring:
                atom_to_aromatic_atoms[idx1].add(idx2)
    
    # 随机选择一个键
    bonds = list(mol.GetBonds())
    selected_bond = random.choice(bonds)
    start_atoms = [selected_bond.GetBeginAtom(), selected_bond.GetEndAtom()]
    
    # 第一步BFS: 拓展距离r内的原子，遇到芳香原子时添加整个芳香环
    visited1 = set()
    dist_dict1 = {}
    queue1 = deque()
    for atom in start_atoms:
        idx = atom.GetIdx()
        visited1.add(idx)
        dist_dict1[idx] = 0
        queue1.append(atom)
    
    while queue1:
        atom = queue1.popleft()
        idx = atom.GetIdx()
        current_dist = dist_dict1[idx]
        if current_dist < r:
            for neighbor in atom.GetNeighbors():
                nidx = neighbor.GetIdx()
                if nidx not in visited1:
                    visited1.add(nidx)
                    dist_dict1[nidx] = current_dist + 1
                    queue1.append(neighbor)
        # 如果当前原子是芳香原子，则添加整个芳香环
        if atom.GetIsAromatic():
            aromatic_atoms = atom_to_aromatic_atoms.get(idx, set())
            for aidx in aromatic_atoms:
                if aidx not in visited1:
                    visited1.add(aidx)
                    # 注意：芳香环原子不加入队列，不分配距离
    
    # 第二步BFS: 搜索距离当前子结构小于d的非芳香原子
    visited2 = set(visited1)
    dist_dict2 = {}
    queue2 = deque()
    for idx in visited1:
        atom = mol.GetAtomWithIdx(idx)
        dist_dict2[idx] = 0
        queue2.append(atom)
    
    while queue2:
        atom = queue2.popleft()
        idx = atom.GetIdx()
        current_dist = dist_dict2[idx]
        for neighbor in atom.GetNeighbors():
            nidx = neighbor.GetIdx()
            if nidx in visited2:
                continue
            if neighbor.GetIsAromatic():
                continue  # 忽略芳香原子
            new_dist = current_dist + 1
            if new_dist < d:  # 距离小于d才添加
                visited2.add(nidx)
                dist_dict2[nidx] = new_dist
                queue2.append(neighbor)
    
    # 生成子结构的SMILES
    smiles = Chem.MolFragmentToSmiles(mol, list(visited2))
    return smiles
import time
if __name__ == "__main__":

    # 重新设置环境变量，指定使用的GPU设备

    predictor_model, graph_functions, template_dicts, template_infos = init_LocalTransform(cfg.predictor)
    def predict_products(reactants):
        if isinstance(reactants, str):
            reactants = [reactants]
        return predict_product_batch(cfg.predictor, reactants, predictor_model, graph_functions, template_dicts, template_infos, verbose = False, sep = False)

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_file",
        type=str,
        default="/root/retro_synthesis/root_align_retro/dataset/filter_weight/uspto_full_test-filter_weight600_rand5_smi.txt",
        help="Path to the input txt file containing the generated molecules.",
    )
    parser.add_argument(
        "--model_ckpt",
        type=str,
        required=True,
        help="Type of used model ckpt.",
    )
    parser.add_argument(
        "--expand",
        action="store_true",
        help="Whether to expand the input molecules by sub_structure.",
    )
    args = parser.parse_args()
    R = 5
    D = 2
    output_prefix1 = f"_{args.model_ckpt}_{f'expand_R{R}_D{D}' if args.expand else 'noexpand'}_retro.txt"
    output_dir_file1 = args.input_file.replace('.txt', f'{output_prefix1}')
    output_prefix2 = f"_{args.model_ckpt}_{f'expand_R{R}_D{D}' if args.expand else 'noexpand'}_retro_forward.txt"
    output_dir_file2 = args.input_file.replace('.txt', f'{output_prefix2}')
    if args.model_ckpt == 'PtoR_FULL':
        model_path = "/root/retro_synthesis/root_align_retro/models/PtoR/USPTO_full_PtoR.pt"   # <<< 换成你的模型路径
    elif args.model_ckpt == 'PtoR_MIT':
        model_path = "/root/retro_synthesis/root_align_retro/models/PtoR/USPTO_MIT_PtoR.pt"
    elif args.model_ckpt == 'PtoR_50K':
        model_path = "/root/retro_synthesis/root_align_retro/models/PtoR/USPTO_50K_PtoR.pt"
    # elif args.model_ckpt == "TMPtoR_FULL":
    #     model_path = "/root/retro_synthesis/root_align_retro/models/USPTO_MIT_TMPtoR_fromUSPTO-FULL_TMP2R_aug10/finetune_average_model_26-30.pt"
    elif args.model_ckpt == "TMPtoR_MIT":
        model_path = "/root/retro_synthesis/root_align_retro/models/USPTO_MIT_TMPtoR_fromUSPTO-MIT_TMP2R_aug10/finetune_average_model_26-30.pt"
    elif args.model_ckpt == "TMPtoR_50K":
        model_path = "/root/retro_synthesis/root_align_retro/models/TMPtoR_fromUSPTO50K_P2R_aug20_smt1/finetune_average_model_26-30.pt"
    elif args.model_ckpt == "TMPtoR_FULL_enhance":
        model_path = "/root/retro_synthesis/root_align_retro/models/TMPtoR_jump3_fromUSPTO_FULL_P2R_aug10_expand-ring/finetune_average_model_26-30.pt"
    elif args.model_ckpt == "TMPtoR_MIT_enhance":
        model_path = "/root/retro_synthesis/root_align_retro/models/TMPtoR_jump3_fromUSPTO_MIT_P2R_aug5_expand-ring/finetune_average_model_26-30.pt"
    elif args.model_ckpt == "TMPtoR_50K_enhance":
        model_path = "/root/retro_synthesis/root_align_retro/models/TMPtoR_jump3_fromUSPTO50K_P2R_aug20_expand-ring/finetune_average_model_26-30.pt"
    elif args.model_ckpt == "TMPtoR_50K_enhance5":
        model_path = "/root/retro_synthesis/root_align_retro/models/TMPtoR_jump5_fromUSPTO50K_P2R_aug20_expand-ring/finetune_average_model_26-30.pt"
    else:
        raise ValueError("Invalid model_ckpt type. Choose from ['PtoR_FULL', 'PtoR_MIT', 'PtoR_50K', 'TMPtoR_FULL', 'TMPtoR_MIT', 'TMPtoR_50K', 'TMPtoR_FULL_enhance', 'TMPtoR_MIT_enhance', 'TMPtoR_50K_enhance']")
    retro = RetroSynthesisDemo(model_path=model_path, gpu_device=0)
    forward = ForwardSynthesisDemo(model_path="/root/retro_synthesis/root_align_retro/models/RtoP/USPTO-MIT_RtoP_mixed.pt", gpu_device=0)
    print(f"Using model checkpoint: {model_path}")
    # Load generated molecules from txt file
    with open(args.input_file, "r") as f:
        data = [f.strip() for f in f.readlines()]
    if args.expand:
        data = [random_substructure(smi, r=R, d=D) for smi in data]
    print(f"Loaded {len(data)} generated molecules from {args.input_file}")
    begin = time.time()
    predict_reactant = retro.inference(data)
    end = time.time()
    print(f"Retro-synthesis inference time: {end-begin:.2f} seconds for {len(data)} molecules")
    print(f"Retro-synthesis done, saving to {output_dir_file1}")

    with open(output_dir_file1, "w") as f:
        for smi in predict_reactant:
            f.write(smi + "\n")
    # Forward synthesis
    begin = time.time()
    predict_product = forward.inference(predict_reactant)
    end = time.time()
    print(f"Forward-synthesis inference time: {end-begin:.2f} seconds for {len(predict_reactant)} molecules")
    print(f"Forward-synthesis done, saving to {output_dir_file2}")
    with open(output_dir_file2, "w") as f:
        for smi in predict_product:
            f.write(smi + "\n")

from rdkit import Chem
import os
import argparse
from tqdm import tqdm
import multiprocessing
import pandas as pd
from rdkit import RDLogger
import re

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)


def smi_tokenizer(smi):
    pattern = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    assert smi == ''.join(tokens)
    return ' '.join(tokens)


def canonicalize_smiles_clear_map(smiles,return_max_frag=True):
    mol = Chem.MolFromSmiles(smiles,sanitize=not opt.synthon)
    if mol is not None:
        [atom.ClearProp('molAtomMapNumber') for atom in mol.GetAtoms() if atom.HasProp('molAtomMapNumber')]
        try:
            smi = Chem.MolToSmiles(mol, isomericSmiles=True)
        except:
            if return_max_frag:
                return '',''
            else:
                return ''
        if return_max_frag:
            sub_smi = smi.split(".")
            sub_mol = [Chem.MolFromSmiles(smiles,sanitize=not opt.synthon) for smiles in sub_smi]
            sub_mol_size = [(sub_smi[i], len(m.GetAtoms())) for i, m in enumerate(sub_mol) if m is not None]
            if len(sub_mol_size) > 0:
                return smi, canonicalize_smiles_clear_map(sorted(sub_mol_size,key=lambda x:x[1],reverse=True)[0][0],return_max_frag=False)
            else:
                return smi, ''
        else:
            return smi
    else:
        if return_max_frag:
            return '',''
        else:
            return ''


def compute_rank(prediction,raw=False,alpha=1.0):
    valid_score = [[k for k in range(len(prediction[j]))] for j in range(len(prediction))]
    invalid_rates = [0 for k in range(len(prediction[0]))]
    rank = {}
    max_frag_rank = {}
    highest = {}
    if raw:
        # no test augmentation
        assert len(prediction) == 1
        for j in range(len(prediction)):
            for k in range(len(prediction[j])):
                if prediction[j][k][0] == "":
                    invalid_rates[k] += 1
            # error detection
            prediction[j] = [i for i in prediction[j] if i[0] != ""]
            for k, data in enumerate(prediction[j]):
                rank[data] = 1 / (alpha * k + 1)
    else:

        for j in range(len(prediction)):
            for k in range(len(prediction[j])):
                # predictions[i][j][k] = canonicalize_smiles_clear_map(predictions[i][j][k])
                if prediction[j][k][0] == "":
                    valid_score[j][k] = opt.beam_size + 1
                    invalid_rates[k] += 1
            # error detection and deduplication
            de_error = [i[0] for i in sorted(list(zip(prediction[j], valid_score[j])), key=lambda x: x[1]) if i[0][0] != ""]
            prediction[j] = list(set(de_error))
            prediction[j].sort(key=de_error.index)
            for k, data in enumerate(prediction[j]):
                if data in rank:
                    rank[data] += 1 / (alpha * k + 1)
                else:
                    rank[data] = 1 / (alpha * k + 1)
                if data in highest:
                    highest[data] = min(k,highest[data])
                else:
                    highest[data] = k
        for key in rank.keys():
            rank[key] += highest[key] * -1e8
    return rank,invalid_rates

import numpy as np
import matplotlib.pyplot as plt
def visualize_results(predictions,opt, accuracy, max_frag_accuracy, invalid_rates, sorted_invalid_rates):
    # Create figure and subplots
    fig, axs = plt.subplots(2, 2, figsize=(18, 12))
    fig.suptitle('Model Performance Metrics: {}'.format(opt.save_file.split('/')[-2] + "--" + opt.save_file.split('/')[-1]) if opt.save_file else 'Model Performance Metrics',
                  fontsize=16, y=1.02)
    
    # Prepare x-axis (top-n values)
    top_n = np.arange(1, opt.n_best + 1)
    display_top_n = [x for x in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 19, 49] if x <= opt.n_best]
    
    # Plot 1: Top-n Accuracy
    axs[0, 0].plot(display_top_n, [accuracy[x-1]/len(predictions)*100 for x in display_top_n], 
                'b-o', linewidth=2, markersize=8)
    axs[0, 0].set_title('Top-n Accuracy', pad=15)
    axs[0, 0].set_xlabel('Top-n', labelpad=10)
    axs[0, 0].set_ylabel('Accuracy (%)', labelpad=10)
    axs[0, 0].set_xlim(1, 10)
    # axs[0, 0].set_ylim(30, 100)
    # y轴最大值是100%，最小值是 min( 60, min(accuracy) / len(predictions) * 100)
    axs[0, 0].set_ylim(min(30, min(accuracy) / len(predictions) * 100), 100)

    axs[0, 0].set_xticks(range(1, 11))
    axs[0, 0].grid(True, alpha=0.3)
    
    # Plot 2: Max Fragment Accuracy
    axs[0, 1].plot(display_top_n, [max_frag_accuracy[x-1]/len(predictions)*100 for x in display_top_n], 
                'g-o', linewidth=2, markersize=8)
    axs[0, 1].set_title('Max Fragment Accuracy', pad=15)
    axs[0, 1].set_xlabel('Top-n', labelpad=10)
    axs[0, 1].set_ylabel('Accuracy (%)', labelpad=10)
    axs[0, 1].set_xlim(1, 10)
    # axs[0, 1].set_ylim(30, 100)
    axs[0, 1].set_ylim(min(30, min(max_frag_accuracy) / len(predictions) * 100), 100)
    axs[0, 1].set_xticks(range(1, 11))
    axs[0, 1].grid(True, alpha=0.3)
    
    # Plot 3: Accuracy Comparison
    axs[1, 1].plot(display_top_n, [accuracy[x-1]/len(predictions)*100 for x in display_top_n], 
                'b-o', label='Overall Accuracy', linewidth=2, markersize=8)
    axs[1, 1].plot(display_top_n, [max_frag_accuracy[x-1]/len(predictions)*100 for x in display_top_n], 
                'g-o', label='Max Frag Accuracy', linewidth=2, markersize=8)
    axs[1, 1].set_title('Accuracy Comparison', pad=15)
    axs[1, 1].set_xlabel('Top-n', labelpad=10)
    axs[1, 1].set_ylabel('Accuracy (%)', labelpad=10)
    axs[1, 1].set_xlim(1, 10)
    # axs[1, 1].set_ylim(30, 100)
    axs[1, 1].set_ylim(min(30, min(accuracy + max_frag_accuracy) / len(predictions) * 100), 100)
    axs[1, 1].set_xticks(range(1, 11))
    axs[1, 1].legend(loc='lower right')
    axs[1, 1].grid(True, alpha=0.3)
    
    # Plot 4: Invalid SMILES Comparison
    axs[1, 0].plot(display_top_n, [invalid_rates[x-1]/len(predictions)/opt.augmentation*100 for x in display_top_n], 
                'r-o', label='Original', linewidth=2, markersize=8)
    axs[1, 0].plot(display_top_n, [sorted_invalid_rates[x-1]/len(predictions)/opt.augmentation*100 for x in display_top_n], 
                'm-o', label='Sorted', linewidth=2, markersize=8)
    axs[1, 0].set_title('Invalid SMILES Comparison', pad=15)
    axs[1, 0].set_xlabel('Top-n', labelpad=10)
    axs[1, 0].set_ylabel('Rate (%)', labelpad=10)
    axs[1, 0].set_xlim(1, 10)
    axs[1, 0].set_xticks(range(1, 11))
    axs[1, 0].legend()
    axs[1, 0].grid(True, alpha=0.3)
    
    # Plot 5: Detailed Accuracy (if available)
    if opt.detailed:
    #     axs[1, 1].plot(display_top_n, [topn_accuracy_chirality[x-1]/total_chirality*100 for x in display_top_n], 
    #                 'c-o', label='With Chirality', linewidth=2, markersize=8)
    #     axs[1, 1].plot(display_top_n, [topn_accuracy_wochirality[x-1]/(len(predictions)-total_chirality)*100 for x in display_top_n], 
    #                 'y-o', label='Without Chirality', linewidth=2, markersize=8)
    #     axs[1, 1].set_title('Chirality Impact', pad=15)
    #     axs[1, 1].set_xlabel('Top-n', labelpad=10)
    #     axs[1, 1].set_ylabel('Accuracy (%)', labelpad=10)
    #     axs[1, 1].set_xlim(1, 10)
    #     axs[1, 1].set_xticks(range(1, 11))
    #     axs[1, 1].legend()
    #     axs[1, 1].grid(True, alpha=0.3)
    # else:
    #     axs[1, 1].axis('off')
        raise NotImplementedError("Detailed accuracy visualization is not implemented yet.")
    
    # Hide the last subplot (6th position)
    
    # Adjust layout to prevent overlap
    plt.tight_layout(pad=3.0)
    
    # Save the figure
    if opt.save_file:
        base_path = os.path.splitext(opt.save_file)[0]
        png_path = f"{base_path}.png"
        plt.savefig(png_path, dpi=300, bbox_inches='tight')
        print(f"Visualization saved to {png_path}")

        csv_path = f"{base_path}.csv"
        top_n = np.arange(1, opt.n_best + 1)
        display_top_n = [x for x in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 19, 49] if x <= opt.n_best]
        data = {
            'Top-n': display_top_n,
            'Accuracy (%)': [accuracy[x-1]/len(predictions)*100 for x in display_top_n],
            'Max Frag Accuracy (%)': [max_frag_accuracy[x-1]/len(predictions)*100 for x in display_top_n],
            'Invalid Rate (%)': [invalid_rates[x-1]/len(predictions)/opt.augmentation*100 for x in display_top_n],
            'Sorted Invalid Rate (%)': [sorted_invalid_rates[x-1]/len(predictions)/opt.augmentation*100 for x in display_top_n]
        }
        df = pd.DataFrame(data)
        df.to_csv(csv_path, index=False)
        print(f"Data saved to {csv_path}")
    else:
        plt.show()
    # Visualize the results
def main(opt):
    print('Reading predictions from file ...')
    with open(opt.predictions, 'r') as f:
        lines = [''.join(line.strip().split(' ')) for line in f.readlines()]
        # lines = lines[:500000]
        print(len(lines))
        data_size = len(lines) // (opt.augmentation * opt.beam_size) if opt.length == -1 else opt.length
        lines = lines[:data_size * (opt.augmentation * opt.beam_size)]
        print("Canonicalizing predictions using Process Number ",opt.process_number)
        pool = multiprocessing.Pool(processes=opt.process_number)
        raw_predictions = pool.map(func=canonicalize_smiles_clear_map,iterable=lines)
        pool.close()
        pool.join()

        predictions = [[[] for j in range(opt.augmentation)] for i in range(data_size)]  # data_len x augmentation x beam_size
        for i, line in enumerate(raw_predictions):
            predictions[i // (opt.beam_size * opt.augmentation)][i % (opt.beam_size * opt.augmentation) // opt.beam_size].append(line)

    print("data size ",data_size)
    print('Reading targets from file ...')
    with open(opt.targets, 'r') as f:
        lines = f.readlines()
        # lines = [''.join(line.strip().split(' ')) for line in f.readlines()]
        print("Origin File Length", len(lines))
        targets = [''.join(lines[i].strip().split(' ')) for i in tqdm(range(0,data_size * opt.augmentation,opt.augmentation))]
        pool = multiprocessing.Pool(processes=opt.process_number)
        targets = pool.map(func=canonicalize_smiles_clear_map, iterable=targets)
        pool.close()
        pool.join()
    ground_truth = targets
    print("Origin Target Lentgh, ", len(ground_truth))
    print("Cutted Length, ",data_size)
    accuracy = [0 for j in range(opt.n_best)]
    topn_accuracy_chirality = [0 for _ in range(opt.n_best)]
    topn_accuracy_wochirality = [0 for _ in range(opt.n_best)]
    topn_accuracy_ringopening = [0 for _ in range(opt.n_best)]
    topn_accuracy_ringformation = [0 for _ in range(opt.n_best)]
    topn_accuracy_woring = [0 for _ in range(opt.n_best)]
    total_chirality = 0
    total_ringopening = 0
    total_ringformation = 0
    atomsize_topk = []
    accurate_indices = [[] for j in range(opt.n_best)]
    max_frag_accuracy = [0 for j in range(opt.n_best)]
    invalid_rates = [0 for j in range(opt.beam_size)]
    sorted_invalid_rates = [0 for j in range(opt.beam_size)]
    unique_rates = 0
    ranked_results = []
    if opt.detailed:
        if not os.path.exists(opt.sources):
            print("Detailed Mode needs the sources.")
            exit(1)
        with open(opt.sources,"r") as f:
            lines = f.readlines()
            ras_src_smiles = [''.join(lines[i].strip().split(' ')) for i in tqdm(range(0,data_size * opt.augmentation,opt.augmentation))]

    for i in tqdm(range(len(predictions))):
        accurate_flag = False
        if opt.detailed:
            chirality_flag = False
            ringopening_flag = False
            ringformation_flag = False
            pro_mol = Chem.MolFromSmiles(ras_src_smiles[i])
            rea_mol = Chem.MolFromSmiles(ground_truth[i][0])
            pro_ringinfo = pro_mol.GetRingInfo()
            rea_ringinfo = rea_mol.GetRingInfo()
            pro_ringnum = pro_ringinfo.NumRings()
            rea_ringnum = rea_ringinfo.NumRings()
            size = len(rea_mol.GetAtoms()) - len(pro_mol.GetAtoms())
            # if (int(ras_src_smiles[i].count("@") > 0) + int(ground_truth[i][0].count("@") > 0)) == 1:
            if ras_src_smiles[i].count("@") > 0 or ground_truth[i][0].count("@") > 0:
                total_chirality += 1
                chirality_flag = True
            if pro_ringnum < rea_ringnum:
                total_ringopening += 1
                ringopening_flag = True
            if pro_ringnum > rea_ringnum:
                total_ringformation += 1
                ringformation_flag = True

        rank, invalid_rate = compute_rank(predictions[i],raw=opt.raw,alpha=opt.score_alpha)
        for j in range(opt.beam_size):
            invalid_rates[j] += invalid_rate[j]
        rank = list(zip(rank.keys(),rank.values()))
        rank.sort(key=lambda x:x[1],reverse=True)
        rank = rank[:opt.n_best]
        ranked_results.append([item[0][0] for item in rank])

        for j, item in enumerate(rank):
            if item[0][0] == ground_truth[i][0]:
                if not accurate_flag:
                    accurate_flag = True
                    accurate_indices[j].append(i)
                    for k in range(j, opt.n_best):
                        accuracy[k] += 1
                    if opt.detailed:
                        atomsize_topk.append((size,j))
                        if chirality_flag:
                            for k in range(j,opt.n_best):
                                topn_accuracy_chirality[k] += 1
                        else:
                            for k in range(j,opt.n_best):
                                topn_accuracy_wochirality[k] += 1
                        if ringopening_flag:
                            for k in range(j,opt.n_best):
                                topn_accuracy_ringopening[k] += 1
                        if ringformation_flag:
                            for k in range(j,opt.n_best):
                                topn_accuracy_ringformation[k] += 1
                        if not ringopening_flag and not ringformation_flag:
                            for k in range(j,opt.n_best):
                                topn_accuracy_woring[k] += 1

        if opt.detailed and not accurate_flag:
            atomsize_topk.append((size,opt.n_best))
        for j, item in enumerate(rank):
            if item[0][1] == ground_truth[i][1]:
                for k in range(j,opt.n_best):
                    max_frag_accuracy[k] += 1
                break
        for j in range(len(rank),opt.beam_size):
            sorted_invalid_rates[j] += 1
        unique_rates += len(rank)

    for i in range(opt.n_best):
        if i in [0,1,2,3,4,5,6,7,8,9,19,49]:
        # if i in range(10):
            print("Top-{} Acc:{:.3f}%, MaxFrag {:.3f}%,".format(i+1,accuracy[i] / data_size * 100,max_frag_accuracy[i] / data_size * 100),
                  " Invalid SMILES:{:.3f}% Sorted Invalid SMILES:{:.3f}%".format(invalid_rates[i] / data_size / opt.augmentation * 100,sorted_invalid_rates[i] / data_size / opt.augmentation * 100))

    print("Unique Rates:{:.3f}%".format(unique_rates / data_size / opt.beam_size * 100))

    if opt.detailed:
        print_topk = [1,3,5,10]
        save_dict = {}
        atomsize_topk.sort(key=lambda x:x[0])
        differ_now = atomsize_topk[0][0]
        topn_accuracy_bydiffer = [0 for _ in range(opt.n_best)]
        total_bydiffer = 0
        for i,item in enumerate(atomsize_topk):
            if differ_now < 11 and differ_now != item[0]:
                for j in range(opt.n_best):
                    if (j+1) in print_topk:
                        save_dict[f'top-{j+1}_size_{differ_now}'] = topn_accuracy_bydiffer[j] / total_bydiffer * 100
                        print("Top-{} Atom differ size {} Acc:{:.3f}%, Number:{:.3f}%".format(j+1,
                                              differ_now,
                                               topn_accuracy_bydiffer[j] / total_bydiffer * 100,
                                               total_bydiffer/data_size * 100))
                total_bydiffer = 0
                topn_accuracy_bydiffer = [0 for _ in range(opt.n_best)]
                differ_now = item[0]
            for k in range(item[1],opt.n_best):
                topn_accuracy_bydiffer[k] += 1
            total_bydiffer += 1
        for j in range(opt.n_best):
            if (j + 1) in print_topk:
                print("Top-{} Atom differ size {} Acc:{:.3f}%, Number:{:.3f}%".format(j + 1,
                      differ_now,
                      topn_accuracy_bydiffer[j] / total_bydiffer * 100,
                      total_bydiffer / data_size * 100))
                save_dict[f'top-{j+1}_size_{differ_now}'] = topn_accuracy_bydiffer[j] / total_bydiffer * 100

        for i in range(opt.n_best):
            if (i+1) in print_topk:
                if total_chirality > 0:
                    print("Top-{} Accuracy with chirality:{:.3f}%".format(i + 1, topn_accuracy_chirality[i] / total_chirality * 100))
                    save_dict[f'top-{i+1}_chilarity'] = topn_accuracy_chirality[i] / total_chirality * 100
                print("Top-{} Accuracy without chirality:{:.3f}%".format(i + 1, topn_accuracy_wochirality[i] / (data_size - total_chirality) * 100))
                save_dict[f'top-{i+1}_wochilarity'] = topn_accuracy_wochirality[i] / (data_size - total_chirality) * 100
                if total_ringopening > 0:
                    print("Top-{} Accuracy ring Opening:{:.3f}%".format(i + 1, topn_accuracy_ringopening[i] / total_ringopening * 100))
                    save_dict[f'top-{i+1}_ringopening'] = topn_accuracy_ringopening[i] / total_ringopening * 100
                if total_ringformation > 0:
                    print("Top-{} Accuracy ring Formation:{:.3f}%".format(i + 1, topn_accuracy_ringformation[i] / total_ringformation * 100))
                    save_dict[f'top-{i+1}_ringformation'] = topn_accuracy_ringformation[i] / total_ringformation * 100
                print("Top-{} Accuracy without ring:{:.3f}%".format(i + 1, topn_accuracy_woring[i] / (data_size - total_ringopening - total_ringformation) * 100))
                save_dict[f'top-{i+1}_wocring'] = topn_accuracy_woring[i] /  (data_size - total_ringopening - total_ringformation)* 100
        print(total_chirality)
        print(total_ringformation)
        print(total_ringopening)
        # df = pd.DataFrame(list(save_dict.items()))
        df = pd.DataFrame(save_dict,index=[0])
        df.to_csv("detailed_results.csv")
    if opt.save_accurate_indices != "":
        with open(opt.save_accurate_indices, "w") as f:
            total_accurate_indices = []
            for indices in accurate_indices:
                total_accurate_indices.extend(indices)
            total_accurate_indices.sort()

            # for index in total_accurate_indices:
            for index in accurate_indices[0]:
                f.write(str(index))
                f.write("\n")

    if opt.save_file != "":
        with open(opt.save_file,"w") as f:
            for res in ranked_results:
                for smi in res:
                    f.write(smi)
                    f.write("\n")
                for i in range(len(res),opt.n_best):
                    f.write("")
                    f.write("\n")

    visualize_results(predictions,opt, accuracy, max_frag_accuracy, invalid_rates, sorted_invalid_rates)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='score.py',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-beam_size', type=int, default=10,help='Beam size')
    parser.add_argument('-n_best', type=int, default=10,help='n best')
    parser.add_argument('-predictions', type=str, required=True,
                        help="Path to file containing the predictions")
    parser.add_argument('-targets', type=str, required=True, help="Path to file containing targets")
    parser.add_argument('-sources', type=str, default="", help="Path to file containing sources")
    parser.add_argument('-augmentation', type=int, default=20)
    parser.add_argument('-score_alpha', type=float, default=1.0)
    parser.add_argument('-length', type=int, default=-1)
    parser.add_argument('-process_number', type=int, default=multiprocessing.cpu_count())
    parser.add_argument('-synthon', action="store_true", default=False)
    parser.add_argument('-detailed', action="store_true", default=False)
    parser.add_argument('-raw', action="store_true", default=False)
    parser.add_argument('-save_file', type=str,default="")
    parser.add_argument('-save_accurate_indices', type=str,default="")

    opt = parser.parse_args()
    print(opt)
    main(opt)

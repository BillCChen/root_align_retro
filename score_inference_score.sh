
# onmt_translate -config train-from-scratch/PtoR/PtoR-50K-aug20-translate.yml
onmt_translate -config pretrain_finetune/finetune/PtoTMPtoR/TMPtoR-50K-aug20-translate.yml

python score.py \
	-beam_size 10 \
	-n_best 10 \
	-augmentation 20 \
	-targets /root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/test/tmp_reactant.txt \
	-process_number 8 \
	-score_alpha 1 \
	-predictions sample/USPTO_50K_TMPtoR_fromUSPTO50KP2R_aug20/finetune_average_model_18-22_all.txt \
	-save_file   sample/USPTO_50K_TMPtoR_fromUSPTO50KP2R_aug20/finetune_average_model_18-22_all_result.txt \
	-synthon 



# onmt_translate -config train-from-scratch/PtoR/PtoR-50K-aug20-translate.yml
# onmt_translate -config pretrain_finetune/finetune/PtoTMPtoR/TMPtoR-50K-aug20-translate.yml
# 50K
# CUDA_VISIBLE_DEVICES=1 onmt_translate -config pretrain_finetune/finetune/PtoTMPtoR/TMPtoR-50K-aug20-translate.yml
# FULL
CUDA_VISIBLE_DEVICES=1 onmt_translate -config pretrain_finetune/finetune/PtoTMPtoR/TMPtoR-FULL-translate.yml
python score.py \
	-beam_size 10 \
	-n_best 10 \
	-augmentation 10 \
	-targets /root/reaction_data/pretrain_aug/USPTO_FULL_PtoTMPtoR_aug10/test/tmp_reactant_mini.txt \
	-process_number 8 \
	-score_alpha 1 \
	-predictions sample/USPTO_FULL_TMPtoR_fromUSPTOFULLP2R_aug10/finetune_model_56-60.txt \
	-save_file   sample/USPTO_FULL_TMPtoR_fromUSPTOFULLP2R_aug10/finetune_model_56-60_result.txt \
	-synthon 


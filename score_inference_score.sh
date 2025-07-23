
# onmt_translate -config train-from-scratch/PtoR/PtoR-50K-aug20-translate.yml
# onmt_translate -config pretrain_finetune/finetune/PtoTMPtoR/TMPtoR-50K-aug20-translate.yml
# 50K
# CUDA_VISIBLE_DEVICES=1 onmt_translate -config pretrain_finetune/finetune/PtoTMPtoR/TMPtoR-50K-aug20-translate.yml
# FULL
CUDA_VISIBLE_DEVICES=1 onmt_translate -config pretrain_finetune/finetune/PtoTMPtoR/TMPtoR-MIT-translate.yml
python score.py \
	-beam_size 10 \
	-n_best 10 \
	-augmentation 10 \
	-targets /root/reaction_data/pretrain_aug/USPTO_MIT_PtoTMPtoR_aug10/test/tmp_reactant.txt \
	-process_number 8 \
	-score_alpha 1 \
	-predictions sample/USPTO_MIT_PtoTMPtoR_fromUSPTOMITP2R_aug10/finetune_model_86-90.txt \
	-save_file   sample/USPTO_MIT_PtoTMPtoR_fromUSPTOMITP2R_aug10/finetune_model_86-90_result.txt \
	-synthon 


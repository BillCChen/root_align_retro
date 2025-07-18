
# onmt_translate -config train-from-scratch/PtoR/PtoR-50K-aug20-translate.yml
# onmt_translate -config pretrain_finetune/finetune/PtoTMPtoR/TMPtoR-50K-aug20-translate.yml
# 50K
# CUDA_VISIBLE_DEVICES=1 onmt_translate -config pretrain_finetune/finetune/PtoTMPtoR/TMPtoR-50K-aug20-translate.yml
# FULL
CUDA_VISIBLE_DEVICES=1 onmt_translate -config pretrain_finetune/finetune/PtoTMPtoR/TMPtoR-FULL-aug10-translate.yml
python score.py \
	-beam_size 10 \
	-n_best 10 \
	-augmentation 10 \
	-targets /root/reaction_data/pretrain_aug/USPTO_full_PtoTMPtoR_aug10/test/tmp_reactant.txt \
	-process_number 8 \
	-score_alpha 1 \
	-predictions output: sample/USPTO_full_TMPtoR_aug10/model_26-30.txt \
	-save_file   output: sample/USPTO_full_TMPtoR_aug10/model_26-30_result.txt \
	-synthon 


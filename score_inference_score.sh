onmt_translate -config train-from-scratch/PtoR/PtoR-50K-aug20-translate.yml

python score.py \
	-beam_size 10 \
	-n_best 10 \
	-augmentation 20 \
	-targets /root/reaction_data/pretrain_aug/USPTO_50K_PtoR_aug20/test/tgt-test_mini.txt \
	-process_number 8 \
	-score_alpha 1 \
	-predictions sample/retrain_USPTO_50K_PtoR_aug20/USPTO_50K_PtoR_rope_36-40.txt \
	-save_file   sample/retrain_USPTO_50K_PtoR_aug20/USPTO_50K_PtoR_rope_36-40_result.txt \
	-synthon 


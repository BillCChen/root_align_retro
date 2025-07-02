onmt_translate -config pretrain_finetune/finetune/PtoTMPtoR/TMPtoR-50K-aug20-translate.yml

python score.py \
	-beam_size 10 \
	-n_best 10 \
	-augmentation 20 \
	-targets /root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/test/tmp_reactant-mini.txt \
	-predictions /root/retro_synthesis/root_align_retro/sample/USPTO_50K_TMPtoR_aug20/finetune_average_model_USPTO_50K_PtoS.txt \
	-process_number 8 \
	-score_alpha 1 \
	-synthon \
	-save_file /root/retro_synthesis/root_align_retro/sample/USPTO_50K_TMPtoR_aug20/finetune_average_model_USPTO_50K_PtoS-results.txt

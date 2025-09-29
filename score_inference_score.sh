# CUDA_VISIBLE_DEVICES=0 onmt_translate -config /root/retro_synthesis/root_align_retro/pretrain_finetune/finetune/PtoTMPtoR/TMPtoR-50K-aug20-translate_beam_search.yml
python score.py \
	-beam_size 10 \
	-n_best 10 \
	-augmentation 20 \
	-targets /root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/test/tmp_reactant_jump3_aug20.txt \
	-process_number 8 \
	-score_alpha 1 \
	-predictions sample/USPTO_50K_TMPtoR_jump3_from_USPTO50K_P2R_aug20/finetune_average_model_21-23.txt \
	-save_file   sample/USPTO_50K_TMPtoR_jump3_from_USPTO50K_P2R_aug20/finetune_average_model_21-23_result.txt \
	-synthon 


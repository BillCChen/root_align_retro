# CUDA_VISIBLE_DEVICES=0 onmt_translate -config /root/retro_synthesis/root_align_retro/pretrain_finetune/finetune/PtoTMPtoR/TMPtoR-50K-aug20-translate_beam_search.yml
python score.py \
	-beam_size 10 \
	-n_best 10 \
	-augmentation 10 \
	-targets /root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/test/tmp_reactant_jump10_aug10.txt \
	-process_number 8 \
	-score_alpha 1 \
	-predictions sample/radius_search/USPTO_50K_TMPtoR_jump10_from_scratch_aug10/finetune_model_20.txt \
	-save_file   sample/radius_search/USPTO_50K_TMPtoR_jump10_from_scratch_aug10/finetune_model_20_result.txt \
	-synthon 


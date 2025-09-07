CUDA_VISIBLE_DEVICES=0 onmt_translate -config /root/retro_synthesis/root_align_retro/pretrain_finetune/finetune/PtoTMPtoR/TMPtoR-50K-aug20-translate.yml
python score.py \
	-beam_size 10 \
	-n_best 10 \
	-augmentation 5 \
	-targets /root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/test/cls/cls_2/tmp_reactant_jump1.txt \
	-process_number 8 \
	-score_alpha 1 \
	-predictions sample/USPTO_50K_TMPtoR_jump1_fromUSPTO50K_P2R_aug5_Moe1/finetune_average_model_26-30.txt \
	-save_file   sample/USPTO_50K_TMPtoR_jump1_fromUSPTO50K_P2R_aug5_Moe1/finetune_average_model_26-30_result.txt \
	-synthon 


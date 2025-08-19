CUDA_VISIBLE_DEVICES=0 onmt_translate -config /root/retro_synthesis/root_align_retro/train-from-scratch/PtoR/PtoR-50K-aug20-translate.yml
python score.py \
	-beam_size 10 \
	-n_best 10 \
	-augmentation 20 \
	-targets /root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/test/cls/cls_2/tgt-test.txt \
	-process_number 8 \
	-score_alpha 1 \
	-predictions sample/USPTO_50K_PtoR_from-scrarch_onUSPTO_50K_aug20_Moe_1/USPTO_50K_PtoR_16-20_cls_2.txt \
	-save_file   sample/USPTO_50K_PtoR_from-scrarch_onUSPTO_50K_aug20_Moe_1/USPTO_50K_PtoR_16-20_cls_2_result.txt \
	-synthon 


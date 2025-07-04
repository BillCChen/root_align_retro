CUDA_VISIBLE_DEVICES=0,1 onmt_translate -config pretrain_finetune/finetune/PtoTMPtoR/TMPtoR-50K-aug20-translate.yml 

python score.py \
	-beam_size 10 \
	-n_best 10 \
	-augmentation 20 \
	-targets /root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/val/tmp_reactant-mini.txt \
	-process_number 8 \
	-score_alpha 1 \
	-predictions sample/USPTO_50K_TMPtoR_aug20_trn_val_tst/finetune_TMP_model_26-30_TMP_product_val.txt \
	-save_file   sample/USPTO_50K_TMPtoR_aug20_trn_val_tst/finetune_TMP_model_26-30_TMP_product_val_result.txt \
	-synthon 


CUDA_VISIBLE_DEVICES=0,1 onmt_translate -config pretrain_finetune/finetune/PtoTMPtoR/TMPtoR-50K-aug20-translate.yml 

python score.py \
	-beam_size 10 \
	-n_best 10 \
	-augmentation 20 \
	-targets /root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/test/tmp_reactant-mini.txt \
	-predictions sample/USPTO_50K_PtoTMP_aug20_fromUSPTO50KP2R/finetune_average_model_11-15.txt \
	-save_file sample/USPTO_50K_PtoTMP_aug20_fromUSPTO50KP2R/finetune_average_model_11-15-results.txt \
	-process_number 8 \
	-score_alpha 1 \
	-synthon 
	

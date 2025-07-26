onmt_translate -config pretrain_finetune/finetune/PtoTMPtoR/TMPtoR-FULL-translate.yml
python score.py \
	-beam_size 10 \
	-n_best 10 \
	-augmentation 10 \
	-targets /root/reaction_data/pretrain_aug/CHORISO_PtoTMPtoR_aug5/test/tmp_reactant.txt \
	-process_number 8 \
	-score_alpha 1 \
	-predictions sample/USPTO_FULL_TMPtoR_fromUSPTOFULLP2R_aug10/finetune_model_26-30_CHORISO_mini.txt \
	-save_file   sample/USPTO_FULL_TMPtoR_fromUSPTOFULLP2R_aug10/finetune_model_26-30_CHORISO_mini_result.txt \
	-synthon 


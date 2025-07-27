# onmt_translate -config pretrain_finetune/finetune/PtoTMPtoR/TMPtoR-CHORISO-translate_BJMU.yml
python score.py \
	-beam_size 10 \
	-n_best 10 \
	-augmentation 10 \
	-targets dataset/CHORISO_PtoTMPtoR_aug5/test/tmp_reactant_mini.txt \
	-process_number 8 \
	-score_alpha 1 \
	-predictions sample/CHORISO_TMPtoR_fromUSPTOFULLP2R_aug10/finetune_model_46-50_CHORISO_mini.txt \
	-save_file   sample/CHORISO_TMPtoR_fromUSPTOFULLP2R_aug10/finetune_model_46-50_CHORISO_mini_result.txt \
	-synthon 


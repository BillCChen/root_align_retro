# example usage:
python score.py \
	-beam_size	<the beam_size when translating, default is 10> \
	-n_best	<the n_best when translating, default is 10> \
	-augmentation <times of augmentation when making predictions> \
	-predictions <the path of the prediction result> \
	-targets <the path of the augmented target> \
	-process_number <number of process you want to use. Higher is faster.> \
	-score_alpha <weighting the number of occurrences and ranking, default is 1> \
	-save_file <save the final prediction resutls>
	-detailed <if you want to see the detailed accuracy like chirality, set it True> \
	-source <is only needed if the detailed is true. the path of augmented source> \
	-synthon <if you want to calculate the accuracy of synthons, set it True>
# for most use PtoR
python score.py \
	-beam_size 10 \
	-n_best 10 \
	-augmentation 20 \
	-targets /root/reaction_data/pretrain_aug/USPTO_50K_PtoR_aug20/test/tgt-test.txt \
	-predictions sample/USPTO_50K_PtoR_aug20/USPTO_50K_PtoR.txt \
	-process_number 8 \
	-score_alpha 1 \
	-synthon \
	-save_file  sample/USPTO_50K_PtoR_aug20/USPTO_50K_PtoR_result.txt
# for most use PtoTMP
python score.py \
	-beam_size 10 \
	-n_best 10 \
	-augmentation 20 \
	-targets /root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/test/tmp_product-mini.txt \
	-predictions /root/retro_synthesis/root_align_retro/sample/USPTO_50K_PtoTMP_aug20/finetune_average_model_26-30.txt \
	-process_number 8 \
	-score_alpha 1 \
	-synthon \
	-save_file /root/retro_synthesis/root_align_retro/sample/USPTO_50K_PtoTMP_aug20/finetune_average_model_26-30-results.txt
# for most use TMPtoR
python score.py \
	-beam_size 10 \
	-n_best 10 \
	-augmentation 20 \
	-targets /root/reaction_data/pretrain_aug/USPTO_50K_PtoTMPtoR_aug20/test/tmp_reactant-mini.txt \
	-predictions /root/retro_synthesis/root_align_retro/sample/USPTO_50K_TMPtoR_aug20/finetune_average_model_26-30.txt \
	-process_number 8 \
	-score_alpha 1 \
	-synthon \
	-save_file /root/retro_synthesis/root_align_retro/sample/USPTO_50K_TMPtoR_aug20/finetune_average_model_26-30-results.txt

#  detailed ~ source : for detailed scoring for chirality and stereochemistry
python score.py \
	-beam_size 10 \
	-n_best 10 \
	-augmentation 20 \
	-targets /root/reaction_data/pretrain_aug/USPTO_50K_PtoR_aug20/test/tgt-test.txt \
	-predictions /root/retro_synthesis/root_align_retro/sample/USPTO_50K_PtoR_aug20/USPTO_50K_PtoR.txt \
	-process_number 8 \
	-score_alpha 1 \
	-save_file /root/retro_synthesis/root_align_retro/sample/USPTO_50K_PtoR_aug20/USPTO_50K_PtoR_result.txt \
	-detailed \
	-source ./dataset/USPTO_50K_PtoR_aug20/test/src-test.txt



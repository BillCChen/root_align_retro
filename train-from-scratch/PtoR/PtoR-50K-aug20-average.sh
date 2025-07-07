onmt_average_models -output  models/TMPtoR_fromUSPTO50K_P2R_aug20/average_model_18-22.pt \
    -m  models/TMPtoR_fromUSPTO50K_P2R_aug20/finetune_model.product-template_step_180000.pt \
        models/TMPtoR_fromUSPTO50K_P2R_aug20/finetune_model.product-template_step_190000.pt \
        models/TMPtoR_fromUSPTO50K_P2R_aug20/finetune_model.product-template_step_200000.pt \
        models/TMPtoR_fromUSPTO50K_P2R_aug20/finetune_model.product-template_step_210000.pt \
        models/TMPtoR_fromUSPTO50K_P2R_aug20/finetune_model.product-template_step_220000.pt

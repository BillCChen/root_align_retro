onmt_average_models -output  ./models/TMPtoR_fromUSPTO50K_P2R/finetune_average_model_11-15.pt \
    -m  models/TMPtoR_fromUSPTO50K_P2R/finetune_model.product-template_step_110000.pt \
        models/TMPtoR_fromUSPTO50K_P2R/finetune_model.product-template_step_120000.pt \
        models/TMPtoR_fromUSPTO50K_P2R/finetune_model.product-template_step_130000.pt \
        models/TMPtoR_fromUSPTO50K_P2R/finetune_model.product-template_step_140000.pt \
        models/TMPtoR_fromUSPTO50K_P2R/finetune_model.product-template_step_150000.pt

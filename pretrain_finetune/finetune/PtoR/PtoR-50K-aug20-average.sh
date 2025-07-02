onmt_average_models -output  ./models/TMPtoR/finefinetune_average_model_6-10.pt \
    -m  models/TMPtoR/finefinetune_model.product-template_step_60000.pt \
        models/TMPtoR/finefinetune_model.product-template_step_70000.pt \
        models/TMPtoR/finefinetune_model.product-template_step_80000.pt \
        models/TMPtoR/finefinetune_model.product-template_step_90000.pt \
        models/TMPtoR/finefinetune_model.product-template_step_10000.pt

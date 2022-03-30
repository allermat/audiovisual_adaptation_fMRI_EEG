clearvars;
close all;
% Comparison of predictions from the spatial (hemifield) and place code
% models. These figures demonstrate that depending on parameter settings
% these two models can predict the same mean activations (first moment) and
% they always predict the same second moment (representational
% dissimilarity). Note, this function plots all other models assessed in 
% the study as well, here we focused only on the spatial (hemifield) and 
% place code models
% Both models predict linearly increasing mean activations: hemifield model
% ipsi-contra tuned neuron ratio is 30-70%, place code tuning functions are
% non-uniformly distributed across azimuth. 
[indiv_activ_table_lin, mean_activ_table_lin] = ...
    figure_model_distributions('hemiPred','lin','placePred','lin');
% Both models predict constant mean activations: hemifield model
% ipsi-contra tuned neuron ratio is 50-50%, place code tuning functions are
% uniformly distributed across azimuth. 
[indiv_activ_table_const, mean_activ_table_const] = ...
    figure_model_distributions('hemiPred','const','placePred','const');
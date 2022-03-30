clearvars;
close all;
% Spatial and decisional uncertainty model predictions extended with
% decisional choice and response time models. 
[indiv_activ_table,mean_activ_table,R_1st,pval_1st] = ...
    figure_model_distributions('plotCorr',true);
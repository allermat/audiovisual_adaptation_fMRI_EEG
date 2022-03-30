function run_generate_model_RDMs()
% Function to generate model RDMs

% The average amount of shift in the group level PSEs 
% calculated as: (PSE_left - PSE_right)/2
shift = 2.3;
seed = 1234;
% Hemifield code model, pre-test
[raw_pre] = population_rate_code_model('propContra',[0.7,0.7],'plotFigure',false,'seed',seed);
% taking only one hemisphere here
RDM_mdl_hemi_pre = squareform(pdist(raw_pre(1:1e5,:)','euclidean'));

% Hemifield code model, post-test with shift
[raw_left_shift] = population_rate_code_model('propContra',[0.7,0.7],...
                                              'shift',-shift,'plotFigure',false,'seed',seed);
[raw_right_shift] = population_rate_code_model('propContra',[0.7,0.7],...
                                                'shift',shift,'plotFigure',false,'seed',seed);
raw_post = cat(2,raw_left_shift(1:1e5,:),raw_right_shift(1:1e5,:));
RDM_mdl_hemi_post_shift = squareform(pdist(raw_post','euclidean'));

% Hemifield code model, post-test no shift
[raw_noshift] = population_rate_code_model('propContra',[0.7,0.7],'plotFigure',false,'seed',seed);
raw_post = cat(2,raw_noshift(1:1e5,:),raw_noshift(1:1e5,:));
RDM_mdl_hemi_post_noshift = squareform(pdist(raw_post','euclidean'));

% Decisional uncertainty model, pre-test (this model does not incorporate
% hemispheres, so it is fine to use the internally generated RDM)
[~,RDM_mdl_dec_pre] = decision_uncertainty_model('plotFigure',false,'seed',seed);

% Decisional uncertainty model, pre-test
[raw_left_shift] = decision_uncertainty_model('shift',-shift,'plotFigure',false,'seed',seed);
[raw_right_shift] = decision_uncertainty_model('shift',shift,'plotFigure',false,'seed',seed);
raw_post = cat(2,raw_left_shift,raw_right_shift);
RDM_mdl_dec_post_shift = squareform(pdist(raw_post','euclidean'));

% Noise model, pre-test
[raw_pre] = random_code_model('plotFigure',false,'seed',seed);
RDM_mdl_noise_pre = squareform(pdist(raw_pre(1:1e5,:)','euclidean'));

% Noise model, post-test 
[raw_noshift] = random_code_model('plotFigure',false,'seed',seed);
raw_post = cat(2,raw_noshift,raw_noshift);
RDM_mdl_noise_post = squareform(pdist(raw_post','euclidean'));

savePath = fullfile(get_path('project'),'results','data','RDM_mdl.mat');
save(savePath,'RDM_mdl_hemi_pre','RDM_mdl_hemi_post_shift',...
     'RDM_mdl_hemi_post_noshift','RDM_mdl_dec_pre','RDM_mdl_dec_post_shift',...
     'RDM_mdl_noise_pre','RDM_mdl_noise_post');
end

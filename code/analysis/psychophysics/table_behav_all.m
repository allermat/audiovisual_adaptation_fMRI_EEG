function behav_table = table_behav_all(dataPath)

% Psychometric function
PF = @PAL_CumulativeNormal;

% Load input data
load(dataPath);

% Data modality
[pathstr, name, ext] = fileparts(dataPath); 
tmp = strsplit(name, '_');
modality = tmp{1};

% Stimulus locations
StimLevels = [-12 -5 -2 0 2 5 12];

% Initialize table
behav_table = table([repmat({'pse'}, 3, 1); repmat({'hit rate'}, 5, 1); ...
    repmat({'dprime'}, 5, 1)], {'pre'; 'postAV'; 'postVA'; ...
    'AVadapt'; 'VAadapt'; 'pre'; 'postAV'; 'postVA'; ...
    'AVadapt'; 'VAadapt'; 'pre'; 'postAV'; 'postVA'}, ...
    'VariableNames',{'measure','condition'});

% Number of subjects
nsubjects = numel(out);

for ss=1:nsubjects
    % Pretest
    NumPos{ss}(1,:) = sum(out(ss).pre.NumPos, 1);
    OutOfNum{ss}(1,:) = sum(out(ss).pre.OutOfNum, 1);
    
    % Posttest
    if ss == 1 && strcmp(modality, 'EEG')
        % postAV adaptation
        NumPos{ss}(2,:) = out(ss).post.NumPos(1,:);
        OutOfNum{ss}(2,:) = out(ss).post.OutOfNum(1,:);
        % postVA adaptation
        NumPos{ss}(3,:) = out(ss).post.NumPos(2,:);
        OutOfNum{ss}(3,:) = out(ss).post.OutOfNum(2,:);
    else
        % postAV adaptation
        NumPos{ss}(2,:) = sum(out(ss).post.NumPos(1:2,:), 1);
        OutOfNum{ss}(2,:) = sum(out(ss).post.OutOfNum(1:2,:), 1);
        % postVA adaptation
        NumPos{ss}(3,:) = sum(out(ss).post.NumPos(3:4,:), 1);
        OutOfNum{ss}(3,:) = sum(out(ss).post.OutOfNum(3:4,:), 1);
    end

    % Initial fit to get parameter estimates
    for i=1:3
        paramsInit(i,:) = fit_PF(PF, StimLevels, NumPos{ss}(i,:), ...
            OutOfNum{ss}(i,:), [], get_PAL_opt_precise);
    end
    
    % Eta parameter for betabinomial model
    eta = 0.15;

    % Fit fuller model (recalibration model)
    thetasInit = [paramsInit(:,1)' mean(paramsInit(:,2)) paramsInit(1,4)];
    paramsIDmatrix = [1 4 5 5; 2 4 5 5; 3 4 5 5];
    mdl_fuller(ss) = multifit_PF(PF, 'fuller model', thetasInit, ...
        paramsIDmatrix, StimLevels, NumPos{ss}, OutOfNum{ss}, [0 0.1], eta);
end

% Group level PSE
params = cat(3, mdl_fuller.params);
pse = squeeze(params(:,1,:));

% Group level hit rate, d-prime
adapt = cat(1, out.adapt);
pre = cat(1, out.pre);
post = cat(1, out.post);
hit_adapt = cat(1, adapt.hit);
hit_pre = cat(1, pre.hit);
hit_post = cat(1, post.hit);
dprime_adapt = cat(1, adapt.dprime);
dprime_pre = cat(1, pre.dprime);
dprime_post = cat(1, post.dprime);

% Add mean values to table
behav_table.mean = [mean(pse, 2); mean(cat(1, hit_adapt.radapt)); ... 
    mean(cat(1, hit_adapt.ladapt)); mean(hit_pre(:,1)); ...
    mean(cat(1, hit_post.radapt)); mean(cat(1, hit_post.ladapt)); ...
    mean(cat(1, dprime_adapt.radapt)); ... 
    mean(cat(1, dprime_adapt.ladapt)); mean(dprime_pre(:,1)); ...
    mean(cat(1, dprime_post.radapt)); mean(cat(1, dprime_post.ladapt))];
   
% Add sem values to table
behav_table.sem = [std(pse, [], 2); std(cat(1, hit_adapt.radapt)); ...
    std(cat(1, hit_adapt.ladapt)); std(hit_pre(:,1)); ...
    std(cat(1, hit_post.radapt)); std(cat(1, hit_post.ladapt)); ...
    std(cat(1, dprime_adapt.radapt)); ...
    std(cat(1, dprime_adapt.ladapt)); std(dprime_pre(:,1)); ...
    std(cat(1, dprime_post.radapt)); std(cat(1, dprime_post.ladapt))] ...
    / sqrt(nsubjects);

% Add subject values to table
for ss=1:nsubjects
    behav_table.(sprintf('sub%d', ss)) = [pse(:,ss); ...
        hit_adapt(ss).radapt; hit_adapt(ss).ladapt; hit_pre(ss); ...
        hit_post(ss).radapt; hit_post(ss).ladapt; ...
        dprime_adapt(ss).radapt; dprime_adapt(ss).ladapt; ...
        dprime_pre(ss); dprime_post(ss).radapt; dprime_post(ss).ladapt];
end
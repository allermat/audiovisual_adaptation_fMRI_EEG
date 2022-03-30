function [pUncorr,hUncorr,pCorrSample,hCorrSample,obsStat] = bootstrpOneSampleTtest(data,nBoot,varargin)
% Second level bootstrap test based on the t-statistic
% 
% USAGE:
%   [pUncorr,hUncorr,pCorrSample,hCorrSample,obsStat] = ...
%       bootstrpOneSampleTtest(data,nboot,varargin)
%   [pUncorr,hUncorr,pCorrSample,hCorrSample,obsStat] = ...
%       bootstrpOneSampleTtest(data,nboot,'Name',Value)
%
% INPUT:
%   Required:
%       data (matrix): nSamples x nSubjects matrix of input data
%       nboot (scalar): number of bootstrap samples
%   'Name'-Value arguments
%       MCPsol (string): solution to the multiple comparison
%           problem, possible values: 'cluster','max','none', 
%           default: 'cluster'
%       clusterThres (scalar): cluster thresholding p-value, 
%           default: 0.05
%       clusterStat (string): cluster statistic, possible values:
%           'maxsum', 'maxsize'
%       alpha (scalar): alpha leve, default: 0.05
%       H0 (scalar): null hypothesis value, default: 0
%
% OUTPUT:
%   pUncorr (column vector): uncorrected sample p-values
%   hUncorr (column vector): uncorrected sample hypothesis test
%       results
%   pCorrSample (column vector): sample p-values corrected for MCP
%   hCorrSample (column vector): sample hypothesis test results 
%       corrected for MCP
%   obsStat (column vector): observed sample statistics
% 

% Copyright(C) 2017, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
validMCPsol = {'cluster','max','bonf_holm','none'};
validClusterStat = {'maxsum','maxsize'};
addRequired(p,'data',...
            @(x) validateattributes(x,{'numeric'},{'ndims',2}));
addRequired(p,'nBoot',...
            @(x) validateattributes(x,{'numeric'},{'scalar','positive'}));
addParameter(p,'MCPsol','cluster',@(x) any(validatestring(x,validMCPsol)));
addParameter(p,'clusterThres',0.05,...
             @(x) validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',1}));
addParameter(p,'clusterStat','maxsize',@(x) any(validatestring(x,validClusterStat)));
addParameter(p,'alpha',0.05,...
             @(x) validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',1}));
addParameter(p,'H0',0,...
             @(x) validateattributes(x,{'numeric'},{'scalar','finite'}));
addParameter(p,'doubleBoot',false,...
             @(x) validateattributes(x,{'logical'},{'scalar'}));

parse(p,data,nBoot,varargin{:});

data = p.Results.data;
nBoot = p.Results.nBoot;
MCPsol = p.Results.MCPsol;
clusterThres = p.Results.clusterThres;
clusterStat = p.Results.clusterStat;
alpha = p.Results.alpha;
H0 = p.Results.H0;
doubleBoot = p.Results.doubleBoot;

% Number of samples and subjects
[nSamples, nSubjects] = size(data);

    % One-sample t-test
    function st = oneSampleTtest(x,mu) 
    
    if isscalar(mu)
        st = sqrt(size(x,1)).*(mean(x)-mu)./std(x);
    elseif size(mu,2) == size(x,2)
        st = sqrt(size(x,1)).*(mean(x)-mu(1,:))./std(x);
    else
        error('Input size mismatch');
    end
    
    end

% Compute observed statistic
obsStat = oneSampleTtest(data',H0);

% Bootstrap sampling, result is nBoot x nsample
sampling = 'random';
if ~doubleBoot
    if nBoot <= nSubjects^nSubjects
        % Monte-Carlo sampling
        bootStat = bootstrp(nBoot,@oneSampleTtest,data', ...
                               repmat(mean(data,2)',nSubjects,1));
    else
        % Exhaustive sampling
        warning('mvpa:bootstrpOneSampleTtest:exhaustiveSampling',...
                ['The number of possible permutations is less than  ' ...
                 'the required, switching to exhaustive sampling!']);
        bootSam = fullfact(repmat(nSubjects,1,nSubjects));
        nBoot = size(bootSam,1);
        
        bootSam = mat2cell(bootSam,ones(nBoot,1),nSubjects);
        bootStat = cellfun(@(x) oneSampleTtest(data(:,x)',mean(data,2)'),...
                              bootSam,'UniformOutput',false);
        bootStat = cell2mat(bootStat);
    end
    
else
    % error('mvpa:bootstrpOneSampleTtest:notImplementedOption',...
          % 'This option has not yet been implemented. ');

    % if nBoot <= nSubjects^nSubjects
        % Bootstrapping, result is nBoot x nsample
    [~,~,bootStat] = boottest(data','mean',0,2,0.05,nBoot,200);    

    % else
    %     % Exhaustive sampling
    %     warning('mvpa:bootstrpOneSampleTtest:exhaustiveSampling',...
    %             ['The number of possible permutations is less than  ' ...
    %              'the required, switching to exhaustive sampling!'])
    %     bootSam = fullfact(repmat(nSubjects,1,nSubjects))';
    %     nBoot = size(bootSam,2);
    %     bootStat = NaN(nBoot,nSamples);
    %     for i = 1:nSamples
    %         [~,bootStat(:,i)] = ibootci([nBoot,round(nBoot/10)],{@oneSampleTtest,data(i,:)',...
    %                                   repmat(mean(data(i,:),2),nSubjects,1)}, ...
    %                                   'bootsam',bootSam);
    %     end
    % end
end


% switch sampling
% %     case 'unique' % sampling unique bootStat values (wrong approach!!! since bootsam ids should be unique cf. 'exhaustive')
% %         for i=2%1:numel(obsStat)
% %             bootStat{i} = bootstrp(nboot,@oneSampleTtest,data(i,:)',repmat(mean(data(i,:)),size(data,2),1));
% %         end
% %         bootStat = cellfun(@unique, bootStat, 'un', 0);
% %         max_samples = max(cellfun(@numel, bootStat));
% %         for i=1:numel(obsStat)
% %             bootStat{i}(end+1:max_samples) = NaN;
% %         end
% %         bootStat = cat(2, bootStat{:});
% %         
%     case 'exhaustive' % exhaustive sampling only in case of very few subjects
%         bootsam = fullfact(repmat(nSubjects, 1, nSubjects))';
%         bootsam = removeRepetition(bootsam); % remove same value sampling
%         nboot = size(bootsam, 2);
%         mn_sample = mean(data, 2);
%         for i=1:nboot
%             for j=1:nSamples
%                 bootStat(i,j) = oneSampleTtest(data(j,bootsam(:,i))', mn_sample(j));
%             end
%         end
        
%     case 'random' % random sampling
%         [bootStat, bootsam] = bootstrp(nboot,@oneSampleTtest,data',repmat(mean(data,2)',nSubjects,1));
% end

pUncorr = arrayfun(@(x) sum(bootStat(:,x) >= obsStat(x))/nBoot,1:nSamples)';
hUncorr = pUncorr <= alpha;

switch MCPsol
    case 'max'
        maxBoot = max(bootStat,[],2);
        % pCorrOmn = sum(maxBoot >= max(obsStat))/nBoot;
        pCorrSample = arrayfun(@(x) sum(maxBoot >= x)/nBoot,obsStat);   
    case 'cluster'
        df = size(data,2)-1;
        crit = tinv(1-clusterThres,df);
        % Right tailed test only
        obsClust = bwconncomp(obsStat >= crit,4);
        % Pre-allocating array for corrected sample p values
        pCorrSample = ones(size(obsStat));
        if isempty(obsClust.PixelIdxList)
            % If there are no clusters, all p values are 1, return
            % here
            pCorrSample = pCorrSample';
            hCorrSample = pCorrSample <= alpha;
            obsStat = obsStat';
            return;
        end
        bootClust = arrayfun(@(i) bwconncomp(bootStat(i,:) >= crit,4),1:nBoot);
        temp = {bootClust.PixelIdxList};
        if strcmp(clusterStat,'maxsum')
            
        elseif strcmp(clusterStat,'maxsize')
            subFun = @(x) max(cellfun(@numel,x));
            maxSizeBoot = cellfun(subFun,temp,'UniformOutput',false);
            maxSizeBoot(cellfun(@isempty,maxSizeBoot)) = {0};
            maxSizeBoot = sort(cell2mat(maxSizeBoot));
            obsClustId = obsClust.PixelIdxList;
            for i = 1:size(obsClustId,2)
                pCorrSample(obsClustId{i}) = sum(maxSizeBoot >= size(obsClustId{i},1))/nBoot;
            end
            
        end
    case 'bonf_holm' % Bonferroni-Holm correction
        [pCorrSample, hCorrSample] = bonf_holm(pUncorr, alpha);
        obsStat = obsStat';
        return
        
    otherwise
        pCorrSample = pUncorr;
end

pCorrSample = pCorrSample';
hCorrSample = pCorrSample <= alpha;
obsStat = obsStat';

end


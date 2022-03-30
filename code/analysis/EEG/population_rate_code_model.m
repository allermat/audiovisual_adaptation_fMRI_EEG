function [y,y_sampl,RDM_euc_raw,RDM_euc_sampl] = population_rate_code_model(varargin)
% Function for population rate code modelling
% 
% DETAILS: The simulation is based on Salminen et al. 2009 and
%   Trapeau & Schonweisner 2015
% 
% USAGE:
%   population_rate_code_model('Name',value)
% INPUT:
%   'Name',value pair arguments:
%       aloc: vector of auditory locations of interest
%       nNeurons: vector of number of neurons for each hemisphere,
%           default = [1e5,1e5];
%       nSamples: vector of number of random samples taken for each
%           hemisphere. Must be the same for both hemispheres. 
%           default = [25,25];
%       propContra: proportion of contralaterally tuned neurons,
%           default = 0.7
%       sigma: standard deviation of tuning curves in degrees
%           default = 67
%       popCenter: center of the populations of tuning curves in
%           degrees, default = 90
%       span: spread of the tuning curves around the center in
%           degrees, default = 20
%       shift: lateral shift in tuning curves in degrees, default = 0
%       latRatio: ratio of overall activity levels across
%           hemispheres left/right. 
%       plotFigure: whether or not to plot output figure
%       seed: seed for the random number generator
%   
% OUTPUT:
%   y: simulated neuron level responses 
%       (nNeurons x nHemispheres) x nConditions
%   y_sampl: simulated average responses sampled from neuron
%       populations - (nSamples x nHemispheres) x nConditions
%   Optional figure with raw and sampled data images and RDMs
%

% Copyright (C) 2017, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
addParameter(p,'aloc',[-12,-5,-2,0,2,5,12],@(x) validateattributes(x,{'numeric'},...
                                                  {'row','finite'}));
addParameter(p,'nNeurons',[1e5,1e5],@(x) validateattributes(x,{'numeric'},...
                                                  {'row','ncols',2,'positive','integer'}));
addParameter(p,'nSamples',[25,25],@(x) validateattributes(x,{'numeric'},...
                                                  {'row','ncols',2,'positive','integer'}));
addParameter(p,'propContra',[0.7,0.7],@(x) validateattributes(x,{'numeric'},...
                                                  {'row','ncols',2,'>=',0,'<=',1}));
addParameter(p,'sigma',67,@(x) validateattributes(x,{'numeric'}, ...
                                                {'scalar','nonnegative','finite'}));
addParameter(p,'popCenter',90,@(x) validateattributes(x,{'numeric'}, ...
                                                  {'scalar','nonnegative','<=',360}));
addParameter(p,'span',20,@(x) validateattributes(x,{'numeric'}, ...
                                                {'scalar','finite','nonnegative'}));
addParameter(p,'shift',0,@(x) validateattributes(x,{'numeric'}, ...
                                                 {'scalar','finite'}));
addParameter(p,'latRatio',1,@(x) validateattributes(x,{'numeric'},...
                                                  {'scalar','finite','nonnegative'}));
addParameter(p,'seed',[],@(x) validateattributes(x,{'numeric'},...
                                                  {'scalar','finite','nonnegative'}));
addParameter(p,'plotFigure',true,@(x) validateattributes(x,{'logical'},{'scalar'}));

parse(p,varargin{:})

aloc = p.Results.aloc;
nNeurons = p.Results.nNeurons;
nSamples = p.Results.nSamples;
propContra = p.Results.propContra;
sigma = p.Results.sigma;
popCenter = p.Results.popCenter;
span = p.Results.span;
shift = p.Results.shift;
latRatio = p.Results.latRatio;
plotFigure = p.Results.plotFigure;
seed = p.Results.seed;

% Further specific input check
nSamples = unique(nSamples);
if numel(nSamples) > 1
    error('population_rate_code_model:invalidInput',...
          'nSamples must be the same for each heimsphere');
end

% Reseeding the random number generator
if isempty(seed)
    rng('shuffle');
else
    rng(seed);
end

propIpsi = 1-propContra;

% Mean of tuning curves uniformly distributed around +90 and -90 degrees
muContra = rand(floor(nNeurons(1)*propContra(1)),1)*span-(span/2)+popCenter;
muIpsi = rand(ceil(nNeurons(1)*propIpsi(1)),1)*span-(span/2)-popCenter;

% Computing neural responses to the sounds left hemisphere. The
% maximum of each tuning curve is 1. 
% fun = @(m) (aloc+shift)*sign(mean(m));
fun = @(m) normpdf(aloc+shift,m,sigma)/normpdf(m,m,sigma);
y = cat(1,arrayfun(fun,muContra,'UniformOutput',false),...
    arrayfun(fun,muIpsi,'UniformOutput',false));
y = cat(1,y{:})/latRatio;

% Neural responses on the right hemisphere
% New sample of random tuning curve means
muContra = rand(floor(nNeurons(2)*propContra(2)),1)*span-(span/2)-popCenter;
muIpsi = rand(ceil(nNeurons(2)*propIpsi(2)),1)*span-(span/2)+popCenter;
temp = cat(1,arrayfun(fun,muContra,'UniformOutput',false),...
           arrayfun(fun,muIpsi,'UniformOutput',false));
y = cat(1,y,cat(1,temp{:}));

% Masking array for hemispheres
hemispheres = cat(1,ones(nNeurons(1),1),ones(nNeurons(2),1)*2);

% Creating random samples of the neurons
y_sampl = NaN(nSamples,size(y,2),2);
for i = 1:2
    temp = mod(randperm(nNeurons(i)),nSamples);
    temp(temp == 0) = nSamples;
    samples = NaN(size(y,1),1);
    samples(hemispheres == i) = temp;
    % Averaging the neural activity within samples
    temp = arrayfun(@(x) mean(y(samples == x,:)),unique(samples(~isnan(samples))),...
                       'UniformOutput',false);
    y_sampl(:,:,i) = cat(1,temp{:});
end
y_sampl = cat(1,y_sampl(:,:,1),y_sampl(:,:,2));
% Computing RDM for raw data
RDM_euc_raw = squareform(pdist(y','euclidean'));

% Computing RDM for sampled data
RDM_euc_sampl = squareform(pdist(y_sampl','euclidean'));

% Plotting
if plotFigure
    figure;
    set(gcf, 'Position', [100 600 1600 300]);

    subplot(1,4,1);
    imagesc(y); hold on;
    line([0,numel(aloc)+1],[nNeurons(1),nNeurons(1)],'Color','w','LineWidth',2);
    set(gca,'Xtick',1:numel(aloc),'XTickLabels',aloc,'Ytick',...
        [0.5*nNeurons(1),1.5*nNeurons(1)],'YtickLabels',...
        {'left hemisphere','right hemisphere'},'YTickLabelRotation',90);
    title('Raw responses, neurons x locations');

    subplot(1,4,2);
    image(rsa.util.scale01(rsa.util.rankTransform_equalsStayEqual(RDM_euc_raw,1)), ...
          'CDataMapping','scaled');
    set(gca,'CLim',[0 1],'CLimMode','manual');
    set(gca,'Xtick',1:numel(aloc),'XTickLabels',num2cell(aloc),...
            'Ytick',1:numel(aloc),'YTickLabels',num2cell(aloc));
    title('RDM of raw responses (Euclidean)');

    subplot(1,4,3);
    imagesc(y_sampl);
    line([0,numel(aloc)+1],[nSamples(1),nSamples(1)],'Color','w','LineWidth',2);
    set(gca,'Xtick',1:numel(aloc),'XTickLabels',aloc,'Ytick',...
        [0.5*nSamples(1),1.5*nSamples(1)],'YtickLabels',...
        {'left hemisphere','right hemisphere'},'YTickLabelRotation',90);
    title('Population responses, samples x locations');

    subplot(1,4,4);
    image(rsa.util.scale01(rsa.util.rankTransform_equalsStayEqual(RDM_euc_sampl,1)), ...
          'CDataMapping','scaled');
    set(gca,'CLim',[0 1],'CLimMode','manual');
    set(gca,'Xtick',1:numel(aloc),'XTickLabels',num2cell(aloc),...
            'Ytick',1:numel(aloc),'YTickLabels',num2cell(aloc));
    title('RDM of population responses (Euclidean)');
end

end
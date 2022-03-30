function varargout = pcm_fit(Y,M,condVec,partVec,varargin)

p = inputParser;
validRunEffects = {'fixed','random'};
validFitModes = {'group','individual'};

addRequired(p,'Y');
addRequired(p,'M');
addRequired(p,'condVec');
addRequired(p,'partVec');
addParameter(p,'runEffect','fixed',@(x) ismember(x,validRunEffects));
addParameter(p,'fitMode','group',@(x) ismember(x,validFitModes));

% Parsing inputs.
parse(p,Y,M,condVec,partVec,varargin{:});

Y = p.Results.Y;
M = p.Results.M;
condVec = p.Results.condVec;
partVec = p.Results.partVec;
runEffect = p.Results.runEffect;
fitMode = p.Results.fitMode;

if strcmp(fitMode,'group')
    [Tgroup,thetaGr,GpredGr] = pcm_fitModelGroup(Y,M,partVec,condVec,...
        'runEffect',runEffect,'fitScale',1);
    [Tcross,thetaCr,GpredCr] = pcm_fitModelGroupCrossval(Y,M,partVec,condVec,...
        'runEffect',runEffect,'groupFit',thetaGr,'fitScale',1);
    varargout = {Tgroup,Tcross,thetaGr,thetaCr,GpredGr,GpredCr};
    
else
    [Tgroup,thetaGr,GpredGr] = pcm_fitModelIndivid(Y,M,partVec,condVec,'runEffect',runEffect);
    [Tcross,~,thetaCr] = pcm_fitModelIndividCrossval(Y,M,partVec,condVec,'runEffect',runEffect);
    varargout = {Tgroup,Tcross,thetaGr,thetaCr,GpredGr};
end

end
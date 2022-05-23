function res = multifit_PF(PF, modelName, thetasInit, paramsIDmatrix, StimLevels, NumPos, OutOfNum, lapseLimits, eta)
% multifit_PF performs multiple fits with model comparison using the
% betabinomial or binomial model provided by David Meijer's custom code

% Define if we want to use betabinomial model or Palamedes's binomial model
if exist('eta', 'var')
    betabinomial = 1;
else
    betabinomial = 0;
end

% Set PAL options to high precision
PAL_opt = get_PAL_opt_precise;

% Make sure that stimulus levels have the right dimensionality
if size(StimLevels, 1) ~= size(NumPos, 1)
    StimLevels = repmat(StimLevels, size(NumPos, 1), 1);
end

% Fit data in a while loop to improve fminsearch convergence
ntries = 0;
exitflag = 0;
if betabinomial
    thetasInit2 = [thetasInit eta];
else
    thetasInit2 = thetasInit;
end
while (ntries < 10) && (exitflag == 0)
    [thetas, negLL, exitflag, output] = PAL_minimize(@DM_PFML_negLLMultiple_Beta, ...
        thetasInit2, PAL_opt, paramsIDmatrix, StimLevels, NumPos, OutOfNum, ...
        PF, lapseLimits, betabinomial);
    ntries = ntries+1;
    thetasInit2 = thetas+1e-4;
end
if ~exitflag
    warning('MLE fit didn''t converge for %s', modelName);
end

% Save results to output
res.exitflag = exitflag;
res.params = thetas(paramsIDmatrix);
res.LL = -negLL;
if betabinomial
    res.eta = thetas(end);
end

% Compute loglikelihood of saturated model and deviance
if strncmp(modelName, 'bootstrap', 9) || strncmp(modelName, 'fuller', 5)
    LLS = DM_PFML_negLLNonParametricMultiple(NumPos, OutOfNum);
    res.Dev = 2 * (LLS - res.LL);
end

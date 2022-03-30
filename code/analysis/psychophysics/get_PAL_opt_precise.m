function PAL_opt = get_PAL_opt_precise

%Optional arguments
PAL_opt = PAL_minimize('options');   %PAL_minimize search options
PAL_opt.TolFun = 1e-12;     %Increase desired precision on LL
PAL_opt.TolX = 1e-12;       %Increase desired precision on parameters
PAL_opt.MaxFunEvals = 50000; %Allow more function evals
PAL_opt.MaxIter = 50000;     %Allow more iterations
PAL_opt.Display = 'off';    %suppress fminsearch messages
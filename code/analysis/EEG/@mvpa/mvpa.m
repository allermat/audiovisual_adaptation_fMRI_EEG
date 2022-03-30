classdef mvpa

    methods (Static)
        
        % Method for computing generalization performance estimates
        mvparesObj = addGenPerfEstimates(mvparesObj)
        
        % Method for adding neurometric functions
        mvparesObj = addNeurometricFunctions(mvparesObj);
        
        % Method for computing recalibration index
        mvparesObj = addRecalIndex(mvparesObj)
        
        % Method for computing various statistics
        mvparesObj = addStats(mvparesObj,var,varargin)
        
        % Method for appending data to an existing MVPA dataset
        appenddata(I)
        
        % Averages mvpares objects
        mvparesObj = averageMvpaRes(I)
        
        % Method for merging mvpares objects
        mvparesObj = mergeMvpaRes(I)
        
        % Prepares an MVPA dataset and saves it on hard drive
        prepdata(I)
        
        % Method for plotting EEG topographies from the MVPA dataset
        showTopographies(filePath)
        
        % Method for generalizing a support vector machine model on
        % timeseries data
        [mvparesObj,filePath] = svmGeneralize(I)
        
        % Method for training a support vector machine model on timeseries 
        % data
        varargout = svmTrain(I)
        
    end
    
    methods (Static) % (Static, Access = private) 
        
        % Method for computing covariance matrix on EEG data
        sigma = compCovEEG(y,X)
        
        % Method for computing residuals on EEG data based on least squares estimate
        r = compResidualsEEG(y,X)
        
        % Method for computing confidence interval of mean
        [ciLow,ciHigh] = confmean(A,dim,alpha)
        
        % Method for groups the examples for cross-validation
        [groupings,rSeed] = assigngroupings(cvScheme,trExamples_info,misc)
        
        % Mehtod for averaging individual trials
        [info,infoAvg,featAvg,seed] = averagetrials(info,feat,condDef,mode,varargin)
        
        % Method for bootstrap based ttest
        [pUncorr,hUncorr,pCorrSample,hCorrSample,obsStat] = bootstrpOneSampleTtest(data,nboot,varargin)
        
        % Method to compute various performance estimates across
        % training and generalization samples.
        varargout = compPerfEstimates(trueLab,predLab)
        
        % Method for computing statistics for generalization performance
        stats = compGenPerfStats(indivData,time)
        
        % Method for computing statistics for recalibration index
        stats = compRecalIndexStats(indivData,time,varargin)
        
        % Method to compute recalibration index
        [RI,RI_err] = compRecalIndex(predLabelInfo,predLabels)
        
        % Method for performing a full K-fold crossvalidation
        [models,scParam,acc,cvMisc] = dokfoldcv(feats,labs,groupings,cvDetails)
        
        % Method for performing a leave one session out crossvalidation
        [models,scParam,acc,cvMisc] = dolosocv(feats,labs,grouping,cvDetails)
        
        % Extracts the necessary data for MVPA for a given subject.
        [feat,info,misc] = extractdata(dataDir,condDef,idStr,tr_method)
        
        % Extracts the required set of features
        feat = extractfeatures(dataM,isExample,trMethod,trLabel,timePoints,varargin)
        
        % Method for fitting neurometric functions
        NMF = fitNeurometricFunctions(predLabels,infoGenExamples,varargin)
        
        % Auxilliari method for fitting neurometric functions
        % across timepoints
        varargout = fitNMFacrossTime(PF,stimLevels,numPos,outOfNum,varargin)
        
        % Performs grid search with two passes on the given parameters with nested cross-validation.
        [bestparam,details] = gridsearch(TR_labs,TR_feats,svm_type,paramCoarse,nFold,nStepsFine)
        
        % Returns the sample closest to some time point in the specified time point set.
        res = indsamplecustom(t,nsample,fs,starttime)
        
        % Converts the libsvm's svmtrain output model converted to a matrix back to sturcture form.
        mdl = mat2mdl(M)
        
        % Converts the libsvm's svmtrain output model to a matrix.
        M = mdl2mat(mdl)
        
        % Method for noise normalization on EEG data
        y = noiseNormalizeEEG(y,s)
        
        % Method for selecting necessary data for mvpa
        [info,feat,isExampleOrig,trialGrouping] = selectdata(dataM,cond,trMethod,trLabel,timePoints,varargin)
        
        % Scales the features according to the specified method.
        [scaledFeats,varargout] = scalefeatures(feats,method,varargin)
        
        % Selects examples matching the given input condition
        isExample = selectexamples(inputConds,condDef,dataInfo)
        
        % Selects the label for the given input label and condition
        labelOut = selectlabel(inputLabel,inputConds)
        
    end
    
end
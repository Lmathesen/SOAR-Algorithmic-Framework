classdef Problem2d < handle
  %MYPROBLEM Summary of this class goes here

  properties
    ListModelType  = {'MNEK','SK','RBF'}; %1 if CPD rbfs, 0 if PD rbfs (CPD: conditionally PD)
    rawY           = NaN;
    x0             = NaN;
    y0             = NaN;
    regParam       = NaN;
    regInterval    = NaN;
    estNoiseVar    = NaN;
    nSamples       = NaN;
    rep            = NaN;
    listSamplingCR = {'OFFmEI','ONmEI','gSTAR'};
    fHandle        = NaN;
    xOpt           = NaN;
    altxOpt        = NaN;
    fOpt           = NaN;
    name           = NaN;
    max_Iter       = NaN;
    min_rep        = NaN;
    BdgIter        = NaN;
    lb             = NaN;
    ub             = NaN;
    dim            = NaN;
    xTrain         = NaN;
    yTrain         = NaN;
    varF           = NaN;
    sampleSet      = NaN;
    sampleSetValues= NaN;
  end

  methods

    %%constructor for the class
    function problem  = Problem(~)
      %intialize PROBLEM problem.lb = problem.lb; problem.ub = problem.ub;
      problem.ListModelType = 'MNEK';
      problem.listSamplingCR = 'OFFmEI'; % used in RBF interpolation models
      min_rep = 10;
      problem.BdgIter = problem.min_rep;
      problem.max_Iter = 100;
    end

    function y = fun(obj,xIn,varargin)
      nPoints = size(xIn,1);
      tol=1e-3;
      idx =  all(all((xIn<=repmat(obj.ub,nPoints,1)+tol)&(xIn>=repmat(obj.lb,nPoints,1)-tol)));

      y = NaN(nPoints,1);

      % check dimension and feasibility
      if (size(xIn,2)~= obj.dim)
        error('Input argument dimensions must agree!') ;
      elseif ~idx
        xIn
        error('Infeasible input')  ;
      else
        for iPoint = 1:nPoints
          y(iPoint) = feval(obj.fHandle,xIn(iPoint,:));
        end
      end
    end
    function [yOpt xOpt] = gridMin(obj,xTrain)
        for i=1:size(xTrain,1)
            for j=1:size(xTrain,2)
                y = obj.fHandle(xTrain(i,j));
            end
        end
        [yOpt xmin] = min(y);
        xOpt = xTrain(xmin,:);
    end
  end
end
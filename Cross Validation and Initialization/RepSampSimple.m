
function [M,mY] = RepSampSimple(problem,nSamps,rep)

dim = problem.dim;
ub  = problem.ub;
lb  = problem.lb;
    M = lhsdesign(nSamps,dim,'criterion','maximin');
    for i = 1:nSamps
        M(i,:) = M(i,:).*(ub - lb) + lb;
    end

yTrue = problem.fun(M);
[problem.fOpt xTemp] = min(yTrue);
problem.xOpt = M(xTemp,:);
mY    = NaN(nSamps,1);
mY    = yTrue;

end


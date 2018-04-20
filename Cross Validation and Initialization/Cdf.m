function CD=Cdf(problem,nSamps)
dim = problem.dim;
ub  = problem.ub;
lb  = problem.lb;
CD = lhsdesign(nSamps,dim,'criterion','maximin');
for i = 1:nSamps
    CD(i,:) = CD(i,:).*(ub - lb) + lb;
end
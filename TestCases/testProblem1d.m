classdef testProblem1d < Problem
        methods
        function obj    = testProblem1d()

        obj.name        = 'Problem';
        obj.dim         = 1;
        obj.lb          = 0;
        obj.ub          = 1;
        obj.xOpt        = 0.7460;
        obj.fOpt        = -11.4510;
        obj.fHandle     = str2func([obj.name,'F1d']);
    end
  end
end


classdef testProblem1d_grlee < Problem
        methods
        function obj    = testProblem1d_grlee()

        obj.name        = 'grlee';
        obj.dim         = 1;
        obj.lb          = .5;
        obj.ub          = 2.5;
        obj.xOpt        = .6;
        obj.fOpt        = -.8;
        obj.fHandle     = str2func([obj.name,'12']);
    end
  end
end


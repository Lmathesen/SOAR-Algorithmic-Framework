classdef testProblem2d_camel6 < Problem2d
        methods
        function obj    = testProblem2d_camel6()

        obj.name        = 'camel6';
        obj.dim         = 2;
        obj.lb          = [-2, -1];
        obj.ub          = [2, 1];
        obj.xOpt        = [0.0898, -0.7126];
        obj.altxOpt     = [0.0898, 0.7126];
        obj.fOpt        = -1.0316;
        obj.fHandle     = str2func([obj.name]);
    end
  end
end


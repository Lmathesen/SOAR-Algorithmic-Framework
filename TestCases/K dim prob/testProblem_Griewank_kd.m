classdef testProblem_Griewank_kd < Problem2d
        methods
        function obj    = testProblem_Griewank_kd()

        obj.name        = 'griewank';
        obj.dim         = 4;
        obj.lb          = ones(1,obj.dim).*(-5);
        obj.ub          = ones(1,obj.dim).*(5);
        obj.xOpt        = zeros(1,obj.dim);
        obj.altxOpt     = NaN;
        obj.fOpt        = 0;
        obj.fHandle     = str2func([obj.name]);
    end
  end
end


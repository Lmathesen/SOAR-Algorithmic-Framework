classdef testProblem_Ackley_kd < Problem2d
        methods
        function obj    = testProblem_Ackley_kd()

        obj.name        = 'ackley';
        obj.dim         = 5;
        obj.lb          = ones(1,obj.dim).*(-32.768);
        obj.ub          = ones(1,obj.dim).*(32.768);
        obj.xOpt        = zeros(1,obj.dim);
        obj.altxOpt     = NaN;
        obj.fOpt        = 0;
        obj.fHandle     = str2func([obj.name]);
    end
  end
end


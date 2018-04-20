classdef testProblem_Rosen_kd < Problem2d
        methods
        function obj    = testProblem_Rosen_kd()

        obj.name        = 'rosen';
        obj.dim         = 5;
        obj.lb          = ones(1,obj.dim).*(-1);
        obj.ub          = ones(1,obj.dim).*(1);
        obj.xOpt        = ones(1,obj.dim);
        obj.altxOpt     = NaN;
        obj.fOpt        = 0;
        obj.fHandle     = str2func([obj.name]);
    end
  end
end


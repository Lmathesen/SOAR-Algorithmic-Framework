classdef testProblemd_sinus4 < Problem4d
        methods
        function obj    = testProblemd_sinus4()

        obj.name        = 'sinus';
        obj.dim         = 4;
        obj.lb          = [0 0 0 0];
        obj.ub          = [pi() pi() pi() pi()];
        obj.xOpt        = [pi()/2 pi()/2 pi()/2 pi()/2];
        obj.fOpt        = -3.5;
        obj.fHandle     = str2func([obj.name]);
    end
  end
end


classdef MCFOPtest
    properties
        sens
        ncoils
        adjoint
        imsize
        phase
    end
    methods
        function obj= set_sens(obj,phasemap)
            obj.phase=phasemap;
        end
    end
    
end

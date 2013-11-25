classdef (Abstract) chebDiscretization 
    
    properties
        domain = []
        dimension = []
        source = []
    end
        
    properties (Dependent)
        numIntervals
    end

    methods
        function n = get.numIntervals(disc)
            n = length(disc.domain) - 1;
        end        

        function t = isempty(disc)
            t = isempty(disc.source);
        end
        
    end
        
    methods (Abstract)
        values = toValues(disc,f)
        f = toFunction(disc,values)
        A = matrix(disc,varargin)
        b = rhs(disc,f,varargin)
    end
    
end

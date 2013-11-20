classdef (Abstract) blockDiscretization 
    
    properties
        domain = [-1 1]
        dimension = []
        source = []
    end
    
    properties (Dependent)
        numIntervals
    end
    
    methods
        function t = isempty(disc)
            t = isempty(disc.source);
        end
        
        function n = get.numIntervals(disc)
            n = length(disc.domain) - 1;
        end        
    end
        
    methods (Abstract)
        %[x,w] = points(disc)   % appropriate?
        values = toValues(disc,f)
        f = toFunction(disc,values)
        A = discretize(disc,dimension)
    end
    
end

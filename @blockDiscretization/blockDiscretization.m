classdef (Abstract) blockDiscretization 
    
    properties
        dimension = []
        source = []
    end
        
    methods
        function t = isempty(disc)
            t = isempty(disc.source);
        end
        
    end
        
    methods (Abstract)
        %[x,w] = points(disc)   % appropriate?
        values = toValues(disc,f)
        f = toFunction(disc,values)
        A = discretize(disc,dimension)
    end
    
end

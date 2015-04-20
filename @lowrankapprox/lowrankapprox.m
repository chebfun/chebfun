classdef lowrankapprox
    
    properties ( Access = public )
       
        cols 
        
        rows 
        
        pivotValues
        
        pivotLocations 
        
        domain = [-1, 1, -1, 1];
        
    end
    
    
    % ABSTRACT METHODS
    methods ( Access = public, Static = false, Abstract = true )

        % Get method.
        val = get(f, prop);

    end
    
    % CONCRETE METHODS (IMPLEMENTED BY THIS ABSTRACT CLASS.)
    methods ( Access = public, Static = false )
        
        
    end
    
end
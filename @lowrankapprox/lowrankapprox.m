classdef lowrankapprox
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        % COLS: column slices used in low rank representation.
        cols 
        % ROWS: row slices used in low rank representation.
        rows
        % PIVOTVALUES: pivot values used in low rank representation.
        pivotValues
        % PIVOTLOCATIONS: pivot locations used in GE.
        pivotLocations
        % DOMAIN: default is [-1,1] x [-1,1].
        domain = [-1 1 -1 1];
    end
    
    
    % ABSTRACT METHODS
    methods ( Access = public, Static = false, Abstract = true )

        % Get method.
        val = get(f, prop);
        
    end
    
    % CONCRETE METHODS (IMPLEMENTED BY THIS ABSTRACT CLASS.)
    methods ( Access = public, Static = false )
        
        % Extract out the low rank approximation:
        varargout = cdr( f )
        
        % Check if objects have the same domain: 
        out = domainCheck(f, g) 
        
    end
    
end
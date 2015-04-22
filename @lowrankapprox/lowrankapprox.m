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
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HIDDEN METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false, Hidden = true )
        
%         % Check to see if domains are equal.
%         out = domainCheck(f, g)
%         
%         % Scale rows and cols of a CHEBFUN2 so that all pivots are 1
%         F = normalizePivots(F)
%         
%         % Normalize the rows and columns of a CHEBFUN2.
%         F = normalizeRowsAndCols(F, p)
%         
%         % Sample Test in constructor. 
%         pass = sampleTest(f, op, tol, flag)
%         
%         % Is a chebfun2 all positive or negative? 
%         [bol, wzero] = singleSignTest(f) 
%         
         % Get the vertical scale of a Chebfun2.
         vscl = vscale(f) 
    end
    
    % CONCRETE METHODS (IMPLEMENTED BY THIS ABSTRACT CLASS.)
    methods ( Access = public, Static = false )
        
        % Extract out the low rank approximation:
        varargout = cdr( f )
        
        % Check if objects have the same domain: 
        out = domainCheck(f, g) 
        
    end
    
end
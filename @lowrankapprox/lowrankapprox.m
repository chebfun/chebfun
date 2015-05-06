classdef lowrankapprox
%LOWRANKAPPROX   Approximate functions on rectangular domains with low rank approximants.
%
%   Abstract class for approximating smooth 2D functions on rectangular
%   domains using low rank approximations. That is, functions are
%   represented in the form:
%
%              f(x,y)  =   sum_j  d_j c_j(y) r_j(x). 
%
% See also CHEBFUN2.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        
        % LOWRANKAPPROX must be objects in CDR form. The columns, C, should
        % be chebfuns or techs, the rows, R, should be chebfuns or techs,
        % and D is a diagonal matrix. 
        
        % COLS: column slices, C, used in low rank representation.
        cols 
        
        % ROWS: row slices, R, used in low rank representation.
        rows
        
        % PIVOTVALUES: the diagonal entries in D.
        pivotValues
        
        % PIVOTLOCATIONS: pivot locations
        pivotLocations
        
        % DOMAIN: default is [-1,1] x [-1,1].
        domain = [-1 1 -1 1];
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ABSTRACT METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = false, Abstract = true )
        
        % Compose method.
        h = compose(f, op, g);
        
        % Constructor.
        g = constructor(f, op, varargin);
        
        % Get method.
        val = get(f, prop);
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HIDDEN METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false, Hidden = true )
        
        % Check to see if domains are equal.
        out = domainCheck(f, g)
        
        % Scale rows and cols of a CHEBFUN2 so that all pivots are 1
        F = normalizePivots(F)
        
        % Normalize the rows and columns of a CHEBFUN2.
        F = normalizeRowsAndCols(F, p)
        
        % Sample Test in constructor. 
        pass = sampleTest(f, op, tol, flag)
      
        % Is a chebfun2 all positive or negative? 
        [bol, wzero] = singleSignTest(f) 

        % Get the vertical scale of a Chebfun2.
        vscl = vscale(f) 
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONCRETE METHODS (IMPLEMENTED BY THIS ABSTRACT CLASS.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
    end
    
end
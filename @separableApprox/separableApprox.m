classdef separableApprox
%SEPARABLEAPPROX   Approximate functions on logically rectangular domains with low rank approximants.
%
%   Abstract class for approximating smooth 2D functions on logically rectangular
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
        
        % SEPARABLEAPPROX must be objects in CDR form. The columns, C, should
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
        
        % Sample method: samples f on a tensor product grid.  Options
        %   X = SAMPLE(F) returns the matrix of values of F on a tensor
        %   product grid.
        %
        %   [U, D, V] = SAMPLE(F) returns the low rank representation of the
        %   values of F on a tensor product grid. X = U * D * V'.
        %
        %   [U, D, V] = SAMPLE(F,M,N) returns the values of F on a M-by-N
        %   tensor product grid.
        varargout = sample(f, varargin);
        
        % Get bivariate expansion coefficients.  Options
        %  X = COEFFS2( F ) returns matrix of coefficients
        %
        % [C, D, R] = COEFFS2( F ) returns a low rank approximation to
        % bivariat coefficients.
        varargout = coeffs2(f);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HIDDEN METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false, Hidden = true )
        
        % Check to see if domains are equal.
        out = domainCheck(f, g)
        
        % Scale rows and cols of a SEPARABLEAPPROX so that all pivots are 1
        F = normalizePivots(F)
        
        % Normalize the rows and columns of a SEPARABLEAPPROX.
        F = normalizeRowsAndCols(F, p)
        
        % Sample Test in constructor. 
        pass = sampleTest(f, op, tol, flag)
      
        % Is a SEPARABLEAPPROX all positive or negative? 
        [bol, wzero] = singleSignTest(f) 

        % Get the vertical scale of a SEPARABLEAPPROX.
        vscl = vscale(f) 
        
    end
    
end

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
        g = constructor(f, OP, varargin);
        
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
        
        % Absolute value of a LOWRANKAPPROX. (f should have one sign in 
        % its domain)
        f = abs( f )
        
        % Extract out the low rank approximation:
        varargout = cdr( f )
        
        % Cholesky factorization of a LOWRANKAPPROX
        R = chol(f)
        
        % Create LOWRANKAPPROX from real and imaginary parts.
        f = complex(f) 
        
        % Complex conjugate of a LOWRANKAPPROX.
        f = conj(f)
        
        % Basic contour plot of LOWRANKAPPROX
        contour(f)
        
        % Cosine of a LOWRANKAPPROX
        f = cos(f)
        
        % Conjugate transpose of a LOWRANKAPPROX.
        f = ctranspose(f)

        % One-dimensional indefinite integral of a LOWRANKAPPROX.
        f = cumsum(f, dim)

        % Partial derivative of a LOWRANKAPPROX.
        f = diff(f, k, dim)
 
        % Evaluate a LOWRANKAPPROX.
        y = feval(f, x)

        % Flip change of row variables for a LOWRANKAPPROX object.
        f = fliplr(f)
        
        % Flip change of column variables for a LOWRANKAPPROX object.        
        f = flipud(f)
        
        % Imaginary part of a LOWRANKAPPROX.
        f = imag(f)

        % True for an empty LOWRANKAPPROX.
        out = isempty(f)

        % Test if LOWRANKAPPROX objects are equal.
        out = isequal(f, g)

        % True for real LOWRANKAPPROX.
        out = isreal(f)
        
        % True for zero LOWRANKAPPROX objects
        out = iszero(f)
        
        % Approximation rank of a LOWRANKAPPROX.
        len = length(f)

        % Maximum slice of a LOWRANKAPPROX.
        [maxVal, maxPos] = max(f)

        % Minimum slice of a LOWRANKAPPROX.
        [minVal, minPos] = min(f)

        % Minimum and maximum slice of a LOWRANKAPPROX.
        [vals, pos] = minandmax(f)

        % Subtraction of two LOWRANKAPPROX objects.
        f = minus(f, g)

        % Left matrix divide for LOWRANKAPPROX objects.
        X = mldivide(A, B)

        % Right matrix divide for a LOWRANKAPPROX.
        X = mrdivide(B, A)

        % Multiplication of LOWRANKAPPROX objects.
        f = mtimes(f, c)

        % Basic surface plot for LOWRANKAPPROX objects.
        varargout = plot(f, varargin)

        % Addition of two LOWRANKAPPROX objects.
        f = plus(f, g)
        
        % Power function of a LOWRANKAPPROX.
        f = power(f, b)

        % QR factorisation of a LOWRANKAPPROX.
        [Q, R] = qr(f)

        % Right array divide for a LOWRANKAPPROX.
        f = rdivide(f, c, pref)

        % Real part of a LOWRANKAPPROX.
        f = real(f)

        % Restrict a LOWRANKAPPROX to a subinterval.
        f = restrict(f, s)

        % Zero curves of a LOWRANKAPPROX in its domain.
        out = roots(f, varargin)

        % Size of a LOWRANKAPPROX.
        [siz1, siz2] = size(f, varargin)

        % One dimensional definite integral of a LOWRANKAPPROX.
        out = sum(f, dim)

        % LOWRANKAPPROX multiplication.
        f = times(f, g, varargin)
         
        % Transpose for LOWRANKAPPROX objects.
        f = transpose(f)

        % Unary minus of a LOWRANKAPPROX.
        f = uminus(f)

        % Unary plus of a LOWRANKAPPROX.
        f = uplus(f)

    end
    
end
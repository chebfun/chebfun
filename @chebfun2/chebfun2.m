classdef chebfun2
%CHEBFUN2   CHEBFUN2 class for representing functions on [a,b]x[c,d].
% 
%   Class for approximating functions defined on finite rectangles. The 
%   functions should be smooth.
%
% CHEBFUN2(F) constructs a CHEBFUN2 object representing the function F on
% [-1,1]x[-1 1]. F can be a string, e.g., 'sin(x.*y)', a function handle, e.g.,
% @(x,y) x.*y + cos(x), or a matrix of numbers. For the first two, F should in
% most cases be "vectorized" in the sense that it may be evaluated at a matrix
% of points and returns a matrix output. 
%
% CHEBFUN2(F, [A B C D]) specifies a rectangle [A B]x[C D] where the 
% function is defined. A, B, C, D must all be finite.
%
% If F is a matrix, A = (a_ij), the numbers aij are used as function values
% at tensor Chebyshev points of the 2nd kind.
% 
% CHEBFUN2(F, k) returns a rank k approximation to F.
%
% CHEBFUN2(F, [m n]) returns a representation of a bivariate polynomial of 
% degree m in x and n in y. The polynomial is compressed in low rank form 
% and the rank k is still determined adaptively (satisfying k<=min(m,n)+1).
%
% CHEBFUN2(F, k, [A B C D]) or CHEBFUN2(F, [m,n], [A B C D]) is nonadaptive in
% rank or degrees, as above, returning a chebfun2 on the domain [A B]x[C D]. 
% 
% CHEBFUN2(F, 'coeffs') where F is matrix uses the matrix as coefficients in 
% a Chebyshev tensor expansion.
%
% The Chebfun2 software system is based on: 
%
% % A. Townsend and L. N. Trefethen, An extension of Chebfun to two dimensions,
% SIAM J. Sci. Comput., 35 (2013), C495-C518.
%
% See also CHEBFUN, CHEBFUN2V.

% Copyright 2014 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun2 information.

% TODO: Improve documentation of input options.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        % COLS: column slices used in low rank representation.
        cols 
        % ROWS: row slices used in low rank representation.
        rows
        % PIVOTVALUES: pivot values used in low rank representation.
        pivotValues
        % PIVOTLOCATIONS: pivot locations used in GE.
        pivotLocations
        % DOMAIN: rectangle of CHEBFUN2, default is [-1,1] x [-1,1].
        domain = [-1 1 -1 1];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function f = chebfun2(varargin)
            % The main CHEBFUN2 constructor!
            
            % Return an empty CHEBFUN:
            if ( (nargin == 0) || isempty(varargin{1}) )
                return
            end
            
            % Call the constructor, all the work is done here:
            f = constructor(f, varargin{:});
            
        end
        
    end        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
         f = conj(f);
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
    %% PRIVATE METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = private, Static = false )
        
        % The main bulk of the CHEBFUN2 constructor:
        g = constructor(g, op, dom, varargin);
        
    end
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % Convert Chebyshev coefficients to values:
        X = coeffs2vals(U); 
        
        % Convert values to Chebyshev coefficients:
        X = vals2coeffs(U); 
        
        % Padua points to tensor grid:
        [C, V, X, Y] = paduaVals2coeffs( F, dom ); 
        
        % Tensor product of Chebyshev points:
        [xx, yy] = chebpts2(nx, ny, domain, kind);
        
        % Outer-product of two chebfuns:
        F = outerProduct(f, g);   
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Private Static methods implemented by CHEBFUN2 class.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = private, Static = true )
        
    end
    
    
end

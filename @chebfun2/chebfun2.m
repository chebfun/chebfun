classdef chebfun2
% CHEBFUN2 CHEBFUN2 class for constructing functions on [a,b]x[c,d].
% 
%   Class for approximating functions defined on finite rectangles. The 
%   functions should be smooth.
%
% CHEBFUN2(F) constructs a CHEBFUN2 object representing the function F on
% [-1,1]x[-1 1]. F can be a string, e.g., 'sin(x.*y)', a function handle, e.g.,
% @(x,y) x.*y + cos(x), or a matrix of values. For the first two, F should in
% most cases be "vectorized" in the sense that it may be evaluated at a matrix
% of points and returns a matrix output.
%
% If F is a matrix, A = (a_ij), the numbers aij are used as function values
% at tensor Chebyshev points of the 2nd kind. 
%
% CHEBFUN2(F, [A B C D]) specifies a rectangle [A B]x[C D] where the 
% function is defined. A, B, C, D must all be finite.
% 
% CHEBFUN2(F, 'coeffs') where F is matrix uses the matrix as coefficients in 
% a Chebyshev tensor expansion.
%
% The Chebfun2 software system is based on: 
%
% % A. Townsend and L. N. Trefethen, An extension of Chebfun to two dimensions,
% SISC, 35 (2013), C495-C518.
%
% See also CHEBFUN, CHEBFUN2V.

% Copyright 2014 by The University of Oxford and The Chebfun2 Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun2 information.

    
    properties ( Access = public )
        % COLS: column slices used in low rank representation
        cols 
        % ROWS: rows slices used in low rank representation
        rows
        % PIVOTVALUES: pivot values used in low rank representation
        pivotValues
        % PIVOTLOCATIONS: pivot locations used in GE
        pivotLocations
        % DOMAIN: rectangle of CHEBFUN2, default is [-1,1] x [-1,1]
        domain = [-1 1 -1 1];
    end
    
    methods
        
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
    
    % Static methods implemented by CHEBFUN class.
    methods ( Static = true )
        
        X = coeffs2vals( U ); 
        
        X = vals2coeffs( U ); 
        
        [xx, yy] = chebpts2(nx, ny, domain);
        
        % Outer-product of two chebfuns.
        F = outerProduct(f, g);   
        
    end

    % Private methods implemented by CHEBFUN2 class.
    methods ( Access = private )
        
    end
    
    % Static private methods implemented by CHEBFUN2 class.
    methods ( Static = true, Access = private )
        
    end
    
    % Methods implemented by CHEBFUN2 class.
    methods
         f = conj(f);
    end
    
end


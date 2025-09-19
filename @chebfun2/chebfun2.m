classdef chebfun2 < separableApprox
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
% If F is a matrix, F = (f_ij), the numbers f_ij are used as function
% values at tensor Chebyshev points of the 2nd kind. CHEBFUN2(F, 'equi')
% assumes the function values f_ij come from an equispaced tensor grid.
% CHEBFUN2(F, 'equix') assumes the values come from a tensor grid that is
% equispaced in x and 2nd-kind Chebyshev in y. CHEBFUN2(F, 'equiy') assumes
% the values come from a tensor grid that is 2nd-kind Chebyshev in x and
% equispaced in y.
% 
% CHEBFUN2(F, k) returns a rank k approximation to F.
%
% CHEBFUN2(F, [m n]) returns a representation of a bivariate polynomial
% with m coefficients in x and n coefficients in y. The polynomial is
% compressed in low rank form and the rank k is still determined adaptively
% (satisfying k<=min(m,n)+1).
%
% CHEBFUN2(F, k, [A B C D]) or CHEBFUN2(F, [m,n], [A B C D]) is nonadaptive in
% rank or degrees, as above, returning a chebfun2 on the domain [A B]x[C D]. 
% 
% CHEBFUN2(F, 'coeffs') where F is a matrix uses the matrix as coefficients
% in a Chebyshev tensor expansion. CHEBFUN2(F, 'coeffsx') uses an expansion
% of coefficients in x and values in y. CHEBFUN2(F, 'coeffsy') uses an
% expansion of values in x and coefficients in y.
%
% CHEBFUN2(F, 'trig') constructs a CHEBFUN2 object representing a smooth
% and periodic function F on [-1,1]x[-1 1]. The resulting CHEBFUN2 is
% represented using a bivariate Fourier expansion. CHEBFUN2(F, 'trigx')
% constructs a CHEBFUN2 that is periodic in x only, represented as a
% bivariate Fourier-Chebyshev expansion. CHEBFUN2(F, 'trigy') constructs a
% CHEBFUN2 that is periodic in y only, represented as a bivariate
% Chebyshev-Fourier expansion.
%
% CHEBFUN2(F, 'periodic')  is the same as CHEBFUN2(F, 'trig').
% CHEBFUN2(F, 'periodicx') is the same as CHEBFUN2(F, 'trigx').
% CHEBFUN2(F, 'periodicy') is the same as CHEBFUN2(F, 'trigy').
%
% CHEBFUN2(F, pref) constructs a CHEBFUN2 object from F with the options
% determined by the CHEBFUNPREF object pref. CHEBFUN2(F, {prefx, prefy})
% uses the CHEBFUNPREF objects prefx in x and prefy in y, where one of
% prefx and prefy may be [].
%
% The Chebfun2 software system is based on: 
%
% A. Townsend and L. N. Trefethen, An extension of Chebfun to two dimensions,
% SIAM J. Sci. Comput., 35 (2013), C495-C518.
%
% See also CHEBFUN, CHEBFUN2V.

% Copyright 2017 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun2 information.

% TODO: Improve documentation of input options.
    
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
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % Convert Chebyshev coefficients to values:
        varargout = coeffs2vals(U, varargin); 
        
        % Convert values to Chebyshev coefficients:
        varargout = vals2coeffs(U, varargin); 
        
        % Padua points to tensor grid:
        [C, V, X, Y] = paduaVals2coeffs( F, dom ); 
        
        % Tensor product of Chebyshev points:
        [xx, yy] = chebpts2(nx, ny, domain, kind);
        
        % Outer-product of two chebfuns:
        F = outerProduct(f, g);  
        
        % Fast spectrally-accurate Poisson solver:
        u = poisson(f, varargin);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Private Static methods implemented by CHEBFUN2 class.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = private, Static = true )
        
    end
    
    
end

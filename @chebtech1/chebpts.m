function [x, w, v, t] = chebpts(n)
%CHEBPTS   Chebyshev points of 1st kind in [-1, 1].
%   CHEBPTS(N) returns N Chebyshev points of the 1st kind in [-1, 1].
%
%   [X, W] = CHEBPTS(N) returns also a row vector of the weights for
%   Clenshaw-Curtis quadrature (computed using [1,2] ).
%
%   [X, W, V] = CHEBPTS(N) returns, in addition to X and W, the weights V
%   for barycentric polynomial interpolation in the Chebyshev points X.
%
%   [X, W, V, T] = CHEBPTS(N) returns also the angles of X.
%
% See also ANGLES, TRIGPTS, LEGPTS, JACPTS, LAGPTS, and HERMPTS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( n == 0 )     % Special case (no points)
    x = []; 
    w = []; 
    v = [];
    t = [];
    
elseif ( n == 1 ) % Special case (single point)
    x = 0; 
    w = 2; 
    v = 1; 
    t = pi/2;
    
else              % General case
    % Chebyshev points using sine function to preserve symmetry:
    x = sin(pi*((-n+1:2:n-1)/(2*n))).';
    
    % Quadrature weights:            
    if ( nargout > 1 ) 
        w = chebtech1.quadwts(n);
    end
    
    % Barycentric weights:
    if ( nargout > 2 )
        v = chebtech1.barywts(n);
    end
    
    % Angles:
    if ( nargout > 3 )
        t = chebtech1.angles(n);
    end
    
end

end

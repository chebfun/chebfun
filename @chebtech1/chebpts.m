function [x, w, v] = chebpts(n)
%CHEBPTS  Chebyshev points of 1st kind in [-1, 1].
%   CHEBPTS(N) returns N Chebyshev points of the 1st kind in [-1, 1].
%
%   [X, W] = CHEBPTS(N) returns also a row vector of the weights for
%   Clenshaw-Curtis quadrature (computed using [1,2] ).
%
%   [X, W, V] = CHEBPTS(N) returns, in addition to X and W, the weights V
%   for barycentric polynomial interoplation in the Chebyshev points X.
%
%   See also LEGPTS, JACPTS, LAGPTS, and HERMPTS.

%   Copyright 2013 by The University of Oxford and The Chebfun Developers. 
%   See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

%   Reference
%
%   [1] Jean-Paul Berrut & Lloyd N. Trefethen, "Barycentric Lagrange
%   Interpolation", SIAM Review 46 (2004), pp 501-517.

if ( n == 0 )     % Special case (no points!)
    x = []; 
    w = []; 
    v = []; 
    
elseif ( n == 1 ) % Special case (single point)
    x = 0; 
    w = 2; 
    v = 1;     
    
else              % General case
    % Chebyshev points 
    x = cos( (2*(n-1:-1:0)+1)*pi/(2*n) ).';
    
    % Quadrature weights            
    if ( nargout > 1 ) 
        w = chebtech1.quadwts(n);
    end
    
    % Barycentric weights (see reference [1])
    if ( nargout > 2 )
        v = chebtech1.barywts(n);
    end
    
end

end


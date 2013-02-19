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

%   [1] Jorg Waldvogel, "Fast construction of the Fejer and Clenshaw-Curtis
%   quadrature rules", BIT Numerical Mathematics 46 (2006), pp 195-202.
%
%   [2] Jorg Waldvogel, www.chebfun.org/and_beyond/programme/slides/wald.pdf
%
%   [3] Jean-Paul Berrut & Lloyd N. Trefethen, "Barycentric Lagrange
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
    x = cos( (2*(n-1:-1:0)+1)*pi/(2*n) ).';  % (Use of sine enforces symmetry.)
    
    % Quadrature weights            
    if ( nargout > 1 ) 
        w = funcheb1.quadwts(n);
    end
    
    % Barycentric weights (see reference [3])
    if ( nargout > 2 )
        v = funcheb1.barywts(n);
    end
    
end

end


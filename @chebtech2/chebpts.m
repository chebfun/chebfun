function [x, w, v, t] = chebpts(n)
%CHEBPTS   Chebyshev points in [-1, 1].
%   CHEBPTS(N) returns N Chebyshev points of the 2nd kind in [-1,1].
%
%   [X, W] = CHEBPTS(N) returns also a row vector of the weights for
%   Clenshaw-Curtis quadrature (computed using Waldvogel's method: [1], [2]).
%
%   [X, W, V] = CHEBPTS(N) returns, in addition to X and W, the barycentric
%   weights V corresponding to the Chebyshev points X. The barycentric weights
%   are normalised to have infinity norm equal to 1 and a positive first entry.
%
%   [X, W, V, T] = CHEBPTS(N) returns also the angles of X.
%
% See also BARY, QUADWTS, BARYWTS, TRIGPTS, LEGPTS, JACPTS, LAGPTS, and HERMPTS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( n == 0 )     % Special case (no points)
    x = []; 
    w = []; 
    v = []; 
    
elseif ( n == 1 ) % Special case (single point)
    x = 0; 
    w = 2; 
    v = 1;     
    
else              % General case
    % Chebyshev points:
    m = n - 1;
    x = sin(pi*(-m:2:m)/(2*m)).';  % (Use of sine enforces symmetry.)
    
    % Quadrature weights:            
    if ( nargout > 1 ) 
        w = chebtech2.quadwts(n);
    end
    
    % Barycentric weights:            
    if ( nargout > 2 )
        v = chebtech2.barywts(n); 
    end
    
    % Angles:
    if ( nargout > 3 )
        t = chebtech2.angles(n);
    end
    
end

end

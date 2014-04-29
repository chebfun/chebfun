function w = quadwts(n)
%QUADWTS   Quadrature weights for Chebyshev points of 2nd kind.
%   QUADWTS(N) returns the N weights for Clenshaw-Curtis quadrature on 2nd-kind
%   Chebyshev points using the algorithm of Waldvogel.
%
% See also CHEBPTS, BARYWTS.

%   References:
%
%   [1] Joerg Waldvogel, "Fast construction of the Fejer and Clenshaw-Curtis
%   quadrature rules", BIT Numerical Mathematics 46 (2006), pp 195-202.
%
%   [2] Joerg Waldvogel, www.chebfun.org/and_beyond/programme/slides/wald.pdf

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

if ( n == 0 )                      % Special case (no points!)
    w = [];
elseif ( n == 1 )                  % Special case (single point)
    w = 2;
else                               % General case
    n = n-1;
    u0 = 1/(n^2 - 1 + mod(n, 2));  % Boundary weights
    L = 0:n-1;                     % Auxiliary vector 1
    r = 2./(1 - 4*min(L, n-L).^2); % Auxiliary vector 2
    w = [ ifft(r-u0), u0 ];        % Clenshaw-Curtis weights
end

end

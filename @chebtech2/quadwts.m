function w = quadwts(n)
%QUADWTS   Quadrature weights for Chebyshev points of 2nd kind.
%   QUADWTS(N) returns the N weights for Clenshaw-Curtis quadrature on 2nd-kind
%   Chebyshev points using the algorithm of Waldvogel.
%
% See also CHEBPTS, BARYWTS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% References:
%   [1] Joerg Waldvogel, "Fast construction of the Fejer and Clenshaw-Curtis
%        quadrature rules", BIT Numerical Mathematics 46 (2006), pp 195-202.
%   [2] Joerg Waldvogel, www.chebfun.org/and_beyond/programme/slides/wald.pdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%  We use a varient of Waldvogel's algorithm, due to Nick Hale. (See his SIAM
%  2013 talk.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( n == 0 )                      % Special case (no points!)
    w = [];
elseif ( n == 1 )                  % Special case (single point)
    w = 2;
else                               % General case
    c = 2./[1, 1-(2:2:(n-1)).^2];  % Exact integrals of T_k (even)
    c = [c, c(floor(n/2):-1:2)];   % Mirror for DCT via FFT
    w = ifft(c);                   % Interior weights
    w([1,n]) = w(1)/2;             % Boundary weights
end

end

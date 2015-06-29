function w = quadwts(n)
%QUADWTS   Quadrature weights for Chebyshev points of 2nd kind.
%   QUADWTS(N) returns the N weights for Clenshaw-Curtis quadrature on 2nd-kind
%   Chebyshev points.
%
% See also CHEBPTS, BARYWTS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
% We use a variant of Waldvogel's algorithm [1], due to Nick Hale. (See below)
% We note this is similar to Greg Von Winkel's approach, which can be found on
% the MathWorks File Exchange.
%
% Let $f(x) = \sum_{k=0}^nc_kT_k(x)$, then\vspace*{-3pt} }
%   I(f) = v.'*c
% where
%   v = \int_{-1}^1T_k(x)dx = { 2/(1-k^2) : k even
%                             { 0         : k odd
%     = v'*inv(TT)*f(x) where TT_{j,k} = T_k(x_j)
%     = (inv(TT)'*v)'*f(x)
% Therefore
%   I(f) = w.'f(x) => w = inv(TT).'*v;
% Here inv(TT).' = inv(TT) is an inverse discrete cosine transform of Type I.
%
% Furthermore, since odd entries in v are zero, can compute via FFT without
% doubling up from N to 2N (though we still need to double up from N/2 to N to 
% facilitate the use of ifft).
%
% References:
%   [1] Joerg Waldvogel, "Fast construction of the Fejer and Clenshaw-Curtis
%       quadrature rules", BIT Numerical Mathematics 46 (2006), pp 195-202.
%   [2] Greg von Winckel, "Fast Clenshaw-Curtis Quadrature", 
%       http://www.mathworks.com/matlabcentral/fileexchange/6911, (2005)
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

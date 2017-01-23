function w = quadwts(n)
%QUADWTS   Quadrature weights for Chebyshev points of 1st kind.
%   QUADWTS(N) returns the N weights for quadrature on 1st-kind Chebyshev grid.
%
% See also CHEBPTS, BARYWTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
% We use a variant of Waldvogel's algorithm [1], due to Nick Hale. (See below)
% We note this is similar to Greg Von Winkel's approach, which can be found on
% the MathWorks File Exchange.
%
% Let $f(x) = \sum_{k=0}^nc_kT_k(x)$, then\vspace*{-3pt} }
%   I(f) = m.'*c
% where
%   m = \int_{-1}^1T_k(x)dx = { 2/(1-k^2) : k even
%                             { 0         : k odd
%     = m'*inv(TT)*f(x) where TT_{j,k} = T_k(x_j)
%     = (inv(TT)'*m)'*f(x)
% Therefore
%   I(f) = w.'f(x) => w = inv(TT).'*m;
% Here inv(TT).' is the discrete cosine transform of type III.
%
% Furthermore, since odd entries in m are zero, can compute via FFT without
% doubling up from N to 2N.
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
    
    m = 2./[1, 1-(2:2:(n-1)).^2];  % Moments - Exact integrals of T_k (even)
    
    % Mirror the vector for the use of ifft: 
    if ( mod(n, 2) )
        c = [m, -m((n+1)/2:-1:2)]; % n is odd
    else
        c = [m, 0, -m(n/2:-1:2)];  % n is even
    end
    v = exp(1i*(0:n-1)*pi/n);      % weight (rotation) vector
    c = c.*v;                      % Apply the weight vector
    w = real(ifft(c));             % Call ifft
    
end

end

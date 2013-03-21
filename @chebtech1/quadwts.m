function w = quadwts(n)
%QUADWTS   Quadrature weights for Chebyshev points of 1st kind.
%   QUADWTS(N) returns the N weights for quadrature on 1st-kind
%   Chebyshev grid.
%
% See also CHEBPTS, BARYWTS.

%   References:
%
%   [1] Joerg Waldvogel, "Fast construction of the Fejer and Clenshaw-Curtis
%   quadrature rules", BIT Numerical Mathematics 46 (2006), pp 195-202.
%
%   [2] Joerg Waldvogel, www.chebfun.org/and_beyond/programme/slides/wald.pdf

%   Copyright 2013 by The University of Oxford and The Chebfun Developers. 
%   See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Quadrature weights for Chebyshev points of 1st kind (also known as
% Fejer nodes of the 1st kind): cos(k*pi/n), where k = 1/2,3/2,...,n-1/2 (see references 
% [1,2]).

if ( n == 1 )                   % Special case (single point)
    w = 2;
else                            % General case
    L = 0:n-1;
    r = 2./(1-4*min(L,n-L).^2); % auxiliary vector
    
    s1 = sign(n/2-L);           % auxiliary vector
    v1 = s1.*r.*exp(1i*pi/n*L);
    w = real(ifft(v1));         % Chebyshev or Fejer-1 weights
end

end

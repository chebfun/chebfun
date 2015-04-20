function g = cumprod( f, dim )
%CUMPROD  Indefinite product integral of a LOWRANKAPPROX. 
%   G = CUMPROD(F) returns the LOWRANKAPPROX G = exp( cumsum(log(F)) )
% 
%   G = CUMPROD(F, DIM) returns the LOWRANKAPPROX G = exp( cumsum(log(F), DIM) )
%
% See also CUMSUM, SUM, PROD.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Compute directly from the definition: 
if ( nargin == 1 )
    g = exp( cumsum( log( f ) ) ); 
else
    g = exp( cumsum( log( f ), dim ) );
end

end

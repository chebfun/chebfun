function g = cumprod( f, varargin )
%CUMPROD  Indefinite product integral of a chebfun2. 
% 
% G = CUMPROD(F) returns the chebfun G = exp( cumsum(log(F)) )
% 
% G = CUMPROD(F, DIM) returns the chebfun G = exp( cumsum(log(F), DIM) )
%
% See also CUMSUM, SUM, PROD.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Just compute directly from the definition: 
if ( isempty( varargin ) )
    g = exp( cumsum( log( f ) ) ); 
else
    dim = varargin{ 1 }; 
    g = exp( cumsum( log( f ), dim ) );
end

end
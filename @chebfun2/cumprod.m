function g = cumprod( f, varargin )
%CUMPROD  Indefinite product integral of a chebfun2. 
% 
% G = CUMPROD(F) returns the chebfun G = exp( cumsum(log(F)) )
% 
% G = CUMPROD(F,DIM) returns the chebfun G = exp( cumsum(log(F),DIM) )
%
% See also CUMSUM, SUM, PROD.

%% 
% Just compute directly from the definition. 
if ( isempty(varargin) )
    g = exp( cumsum( log( f ) ) ); 
else
    dim = varargin{ 1 }; 
    g = exp( cumsum( log( f ), dim ) );
end

end
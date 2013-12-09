function v = prod(f,varargin)
%PROD  Product integral of a chebfun2. 
% 
%   prod(F) returns the chebfun exp(  sum(log(F)) )
% 
%   prod(F,DIM) returns the chebfun exp( sum(log(F),DIM) )
% 
% See also CUMPROD.

% Just use the definition.

if ( nargin == 1 ) 
    dim == 1; 
else
    dim = varargout{1}; 
end

if ( dim == 1 )  % default to y direction. 
    v = exp( sum( log( f ) ) ); 
elseif ( dim == 2 ) 
    v = exp( sum( log(f), dim ) );
else
    error('CHEBFUN2:PROD:DIM','Unrecognized dimension.');
end

end
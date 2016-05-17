function v = prod(f, varargin)
%PROD     Product integral of a SEPARABLEAPPROX. 
%   PROD(F) returns exp( sum(log(F)) )
% 
%   PROD(F, DIM) returns the chebfun exp( sum(log(F), DIM) )
% 
% See also CUMPROD.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 ) 
    dim = 1; 
else
    dim = varargin{1}; 
end

if ( numel( dim ) ~= 1 )
    error('CHEBFUN:SEPARABLEAPPROX:prod:dim', 'DIM should be either 1 or 2.');
end

v = exp( sum( log(f), dim ) );

end

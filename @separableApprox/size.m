function varargout = size( f , dim) %#ok<INUSL>
% SIZE       Size of a SEPARABLEAPPROX
%   D = SIZE(F) returns the two-element row vector D = [inf,inf].
%
%   [M, N] = SIZE(F) returns M = inf and N = inf.
%
%   M = SIZE(F, DIM) returns the dimension specified by the scalar DIM, which is
%   always inf.

% Copyright 2017 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun2 information.

% NOTE: The size of a SEPARABLEAPPROX object does not depend on f!

if ( ( nargin == 1 ) && ( nargout <= 1 ) )
    varargout = { [Inf, Inf] };
elseif ( ( nargout == 2 ) && ( nargout <=1 ) )
    varargout = { Inf, Inf };
elseif ( ( nargin == 2 ) && ( ( dim == 1 ) || ( dim == 2 ) ) )
    varargout = { Inf };
else
    error('CHEBFUN:SEPARABLEAPPROX:size:outputs', 'Too many output arguments.');
end

end

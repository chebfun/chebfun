function varargout = min(varargin)
%MIN   Minimum value of a CHEBFUN2 in one direction.
%   MIN(f) returns a chebfun representing the minimum of the CHEBFUN2 along the
%   y direction, i.e, MIN(f) = @(x) max( f ( x, : ) )
%
%   MIN(f, [], dim) returns a CHEBFUN representing the minimum of f along the
%   DIM direction. If DIM = 1 is along the y-direction and DIM = 2 is along the
%   x-direction.
%
%   WARNING: This function is not always accurate to full machine precision.
%
%   For the global minimum use MIN2.
%
% See also MAX, MAX2, MIN2, MINANDMAX2.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = min@separableApprox(varargin{:});

end

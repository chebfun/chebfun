function varargout = max(varargin)
%MAX   Maximum value of a CHEBFUN2 in one direction.
%   MAX(f) returns a chebfun representing the maximum of the CHEBFUN2 along the
%   y direction, i.e, MAX(f) = @(x) max( f ( x, : ) )
%
%   MAX(f, [], dim) returns a CHEBFUN representing the maximum of f along the
%   DIM direction. If DIM = 1 is along the y-direction and DIM = 2 is along the
%   x-direction.
%
%   WARNING: This function is not always accurate to the expected precision.
%
%   For the global maximum use MAX2.
%
% See also MIN, MAX2, MIN2, MINANDMAX2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = max@separableApprox(varargin{:});

end

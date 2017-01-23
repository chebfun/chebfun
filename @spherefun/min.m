function varargout = min(varargin)
%MIN   Minimum value of a SPHEREFUN in one direction.
%   MIN(f) returns a chebfun representing the minimum of the SPHEREFUN along the
%   latitude direction, i.e, MIN(f) = @(lambda) max( f ( lambda, : ) )
%
%   MIN(f, [], dim) returns a CHEBFUN representing the minimum of f along the
%   DIM direction. If DIM = 1 is along the latitude-direction and DIM = 2 is along the
%   longitude-direction.
%
%   WARNING: This function is not always accurate to the expected precision.
% 
%   For the global minimum use MIN2.
%
% See also SPHEREFUN/MAX, SPHEREFUN/MAX2, SPHEREFUN/MIN2, SPHEREFUN/MINANDMAX2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = min@separableApprox(varargin{:});

end

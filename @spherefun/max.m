function varargout = max(varargin)
%MAX   Maximum value of a SPHEREFUN in one direction.
%   MAX(f) returns a chebfun representing the maximum of the SPHEREFUN along the
%   latitude direction, i.e, MAX(f) = @(lambda) max( f ( lambda, : ) )
%
%   MAX(f, [], dim) returns a CHEBFUN representing the maximum of f along
%   the DIM direction. If DIM = 1 is along the latitude-direction and DIM =
%   2 is along the longitude-direction.
%
%   WARNING: This function is not always accurate to the expected precision.
% 
%   For the global maximum use MAX2.
%
% See also SPHEREFUN/MIN, SPHEREFUN/MAX2, SPHEREFUN/MIN2, SPHEREFUN/MINANDMAX2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = max@separableApprox(varargin{:});

end

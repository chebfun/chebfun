function varargout = max(varargin)
%MAX   Maximum value of a DISKFUN in one direction.
%   MAX(f) returns a chebfun representing the maximum of the DISKFUN along the
%   angular direction, i.e, MAX(f) = @(theta) max( f ( theta, : ) )
%
%   MAX(f, [], dim) returns a CHEBFUN representing the maximum of f along
%   the DIM direction. If DIM = 1 is along the angular-direction and DIM =
%   2 is along the radial-direction.
%
%   WARNING: This function is not always accurate to the expected precision.
% 
%   For the global maximum use MAX2.
%
% See also DISKFUN/MIN, DISKFUN/MAX2, DISKFUN/MIN2, DISKFUN/MINANDMAX2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = max@separableApprox(varargin{:});

end

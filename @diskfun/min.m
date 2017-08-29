function varargout = min(varargin)
%MIN   Minimum value of a DISKFUN in one direction.
%   MIN(f) returns a chebfun representing the minimum of the DISKFUN along 
%   the latitude direction, i.e, MIN(f) = @(lambda) max( f ( lambda, : ) ).
%
%   MIN(f, [], dim) returns a CHEBFUN representing the minimum of f along
%   direction DIM. DIM = 1 is used to compute max along the radial-direction
%   and DIM = 2 along the angular-direction.
%
%   WARNING: This function is not always accurate to the expected precision.
% 
%   For the global minimum use MIN2.
%
% See also DISKFUN/MAX, DISKFUN/MAX2, DISKFUN/MIN2, DISKFUN/MINANDMAX2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = min@separableApprox(varargin{:});

end
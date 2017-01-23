function varargout = iszero(varargin)
%ISZERO   Check if a DISKFUN is identically zero on its domain.
%   OUT = ISZERO( F ) return 1 if the DISKFUN is exactly the zero function, and
%   0 otherwise. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = iszero@separableApprox(varargin{:});

end

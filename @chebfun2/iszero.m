function varargout = iszero(varargin)
%ISZERO   Check if a CHEBFUN2 is identically zero on its domain.
%   OUT = ISZERO( F ) return 1 if the CHEBFUN2 is exactly the zero function, and
%   0 otherwise.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = iszero@separableApprox(varargin{:});

end

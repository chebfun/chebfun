function varargout = sqrt(varargin)
%SQRT   Square root.
%   SQRT(F) returns the square root of a positive CHEBFUN2 F.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = sqrt@separableApprox(varargin{:});

end

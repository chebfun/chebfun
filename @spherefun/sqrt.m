function varargout = sqrt(varargin)
%SQRT   Square root.
%   SQRT(F) returns the square root of a positive SPHEREFUN F.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = sqrt@separableApprox(varargin{:});

end

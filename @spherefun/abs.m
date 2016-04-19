function varargout = abs(varargin)
%ABS Absolute value of a SPHEREFUN.
%   ABS(F) returns the absolute value of a SPHEREFUN. This function does not work
%   if the function passes through or becomes numerically close to zero.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = abs@separableApprox(varargin{:});

end

function varargout = abs(varargin)
%ABS Absolute value of a SPHEREFUN2.
%   ABS(F) returns the absolute value of a SPHEREFUN2. This function does not work
%   if the function passes through or becomes numerically close to zero.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = abs@separableApprox(varargin{:});
end

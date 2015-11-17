function varargout = min2(varargin)
%MIN2   Global minimum of a CHEBFUN2.
%   Y = MIN2(F) returns the global minimum of F over its domain.
%
%   [Y, X] = MIN2(F) returns the global minimum in Y and its location in X.
%
%  This command may be faster if the OPTIMIZATION TOOLBOX is installed.
%
% See also MAX2, MINANDMAX2.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = min2@separableApprox(varargin{:});

end

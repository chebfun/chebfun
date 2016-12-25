function varargout = max2(varargin)
%MAX2   Global maximum of a CHEBFUN2.
%   Y = MAX2(F) returns the global maximum of F over its domain.
%
%   [Y, X] = MAX2(F) returns the global maximum in Y and its location X.
%
%  This command may be faster if the OPTIMIZATION TOOLBOX is installed.
%
% See also MIN2, MINANDMAX2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = max2@separableApprox(varargin{:});

end

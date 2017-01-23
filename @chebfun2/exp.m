function varargout = exp(varargin)
%EXP  Exponential of a CHEBFUN2
%   EXP(F) returns the exponential of a CHEBFUN2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = exp@separableApprox(varargin{:});

end

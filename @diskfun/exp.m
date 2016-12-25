function varargout = exp(varargin)
%EXP  Exponential of a DISKFUN
%   EXP(F) returns the exponential of a DISKFUN. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = exp@separableApprox(varargin{:});

end

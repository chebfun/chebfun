function varargout = sin(varargin)
%SIN   Sine of a CHEBFUN2.
%
% See also SPHEREFUN/SINH and SPHEREFUN/COS

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = sin@separableApprox(varargin{:});

end

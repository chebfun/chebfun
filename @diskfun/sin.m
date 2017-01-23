function varargout = sin(varargin)
%SIN   Sine of a DISKFUN.
%
% See also DISKFUN/SINH and DISKFUN/COS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = sin@separableApprox(varargin{:});

end
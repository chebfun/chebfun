function varargout = tan(varargin)
%TAN   Tangent of a DISKFUN.
%
% See also DISKFUN/SIN, DISKFUN/COS, DISKFUN/TANH

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = tan@separableApprox(varargin{:});

end

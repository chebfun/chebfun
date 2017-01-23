function varargout = mrdivide(varargin)
% /   Right scalar divide for DISKFUN objects.
%   F/C divides the DISKFUN F by a scalar C.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = mrdivide@separableApprox(varargin{:});

end
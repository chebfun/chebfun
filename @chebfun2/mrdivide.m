function varargout = mrdivide(varargin)
% /   Right scalar divide for CHEBFUN2 objects.
%
%   F/C divides the CHEBFUN2 F by a scalar C.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = mrdivide@separableApprox(varargin{:});

end

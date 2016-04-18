function varargout = mrdivide(varargin)
% /   Right scalar divide for SPHEREFUN2 objects.
%
%   F/C divides the SPHEREFUN2 F by a scalar C.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = mrdivide@separableApprox(varargin{:});
end

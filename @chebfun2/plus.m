function varargout = plus(varargin)
%+   Plus for CHEBFUN2 objects.
%
% F + G adds F and G. F and G can be scalars or CHEBFUN2 objects.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = plus@separableApprox(varargin{:});

end

function varargout = uminus(varargin)
%UMINUS   Unary minus for a CHEBFUN2.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = uminus@separableApprox(varargin{:});

end

function varargout = uminus(varargin)
%UMINUS   Unary minus for a SPHEREFUN. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = uminus@separableApprox(varargin{:});

end

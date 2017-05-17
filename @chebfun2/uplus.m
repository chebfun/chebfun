function varargout = uplus(varargin)
%UPLUS   Unary plus for a CHEBFUN2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = uplus@separableApprox(varargin{:});

end

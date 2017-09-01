function varargout = power(varargin)
%.^	     Pointwise power of a SPHEREFUN. 
%
% F.^G returns a SPHEREFUN F to the scalar power G, a scalar F to the
% SPHEREFUN power G, or a SPHEREFUN F to the SPHEREFUN power G.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = power@separableApprox(varargin{:});

end

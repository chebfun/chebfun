function varargout = power(varargin)
%.^	     Pointwise power of a SPHEREFUN2. 
%
% F.^G returns a SPHEREFUN2 F to the scalar power G, a scalar F to the
% SPHEREFUN2 power G, or a SPHEREFUN2 F to the SPHEREFUN2 power G.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = power@separableApprox(varargin{:});
end

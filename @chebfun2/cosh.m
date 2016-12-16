function varargout = cosh(varargin)
%COSH   Hyperbolic cosine of a CHEBFUN2.
%   COSH(F) returns the hyperbolic cosine of F.
%
% See also SINH, COS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = cosh@separableApprox(varargin{:});

end

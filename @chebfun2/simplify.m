function varargout = simplify(varargin)
% Simplify a CHEBFUN2
%
% F = SIMPLIFY( F ) compressed the representation of F to one that is
% numerically the same, but requires fewer parameters to store. Currently this
% simplifies the polynomial degree of F, but not the rank.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = simplify@separableApprox(varargin{:});

end

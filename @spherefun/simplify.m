function varargout = simplify(varargin)
% Simplify a SPHEREFUN
% 
% F = SIMPLIFY( F ) compressed the representation of F to one that i
% numerically the same, but requires fewer parameters to store. Currently
% this simplifies the polynomial degree of F, but not the rank.
%
% F = SIMPLIFY(F, TOL) does the same as above but uses the scalar TOL instead
% of the default simplification tolerance as the relative threshold level.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = simplify@separableApprox(varargin{:});

end

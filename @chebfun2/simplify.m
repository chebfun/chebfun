function varargout = simplify(varargin)
% Simplify a CHEBFUN2.
%
% F = SIMPLIFY( F ) compresses the representation of F to one that is
% numerically the same, but requires fewer parameters to store. This
% simplifies the bivariate polynomial degree of F, but not its rank.
%
% F = SIMPLIFY(F, TOL) does the same as SIMPLIFY( F ) but uses the scalar 
% TOL instead of the default simplification tolerance as the relative 
% threshold level for compression.
%
% F = SIMPLIFY(F, 'rank') compresses the rank of the representation for F 
% to one that is numerically the same. 
%
% F = SIMPLIFY(F, TOL, 'rank') does the same as SIMPLIFY(F, 'rank') but 
% uses the scalar TOL instead of the default simplification tolerance as 
% the relative threshold level for compression.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = simplify@separableApprox(varargin{:});

end

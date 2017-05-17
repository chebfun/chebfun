function varargout = cdr(varargin)
%CDR decomposition of a CHEBFUN2.
%   [C,D,R] = CDR(F) produces a diagonal matrix D of size length(F) by length(F)
%   and quasimatrices C and R of size inf by length(F) such that f(x,y) = C(y,:)
%   * D * R(x,:)'.
%
%   D = CDR(F) returns a vector containing the pivot values used in the
%   construction of F.
%
% See also PIVOTS, SVD.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = cdr@separableApprox(varargin{:});

end

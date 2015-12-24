function varargout = separableApprox(varargin)
%CHEBFUN2   Approximate functions on logically rectangular domains with low rank approximants.
%
%   Abstract class for approximating smooth 2D functions on logically rectangular
%   domains using low rank approximations. That is, functions are
%   represented in the form:
%
%              f(x,y)  =   sum_j  d_j c_j(y) r_j(x).
%
% See also CHEBFUN2.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = separableApprox@separableApprox(varargin{:});

end

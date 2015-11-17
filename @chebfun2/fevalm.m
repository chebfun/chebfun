function varargout = fevalm(varargin)
% FEVALM   Evaluate a CHEBFUN2.
%
% Z = FEVALM(F, X, Y) returns a matrix of values Z of size length(X)-by-length(Y).
% X and Y should be vectors of doubles. This is equivalent to making a meshgrid
% of the vectors X and Y and then using FEVAL to evaluate at that grid.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = fevalm@separableApprox(varargin{:});

end

function varargout = size(varargin)
% SIZE       Size of a CHEBFUN2
%   D = SIZE(F) returns the two-element row vector D = [inf,inf].
%
%   [M, N] = SIZE(F) returns M = inf and N = inf.
%
%   M = SIZE(F, DIM) returns the dimension specified by the scalar DIM, which is
%   always inf.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = size@separableApprox(varargin{:});

end

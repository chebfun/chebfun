function varargout = diag(varargin)
%DIAG(F)   Diagonal of a SPHEREFUN2.
%   G = DIAG(F) returns the CHEBFUN representing g(x) = f(x, x).
%
%   G = diag(F,C) returns the CHEBFUN representing g(x) = f(x, x+c).
%
% See also TRACE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = diag@separableApprox(varargin{:});
end

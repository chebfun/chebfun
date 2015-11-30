function varargout = minandmax2(varargin)
%MINANDMAX2   Find global minimum and maximum of a CHEBFUN2.
%   Y = minandmax2(F) returns the minimum and maximum value of a CHEBFUN2 over
%   its domain. Y is a vector of length 2 such that Y(1) = min(f(x,y)) and Y(2)
%   = max(f(x,y)).
%
%   [Y, X] = minandmax2(F) also returns the position of the minimum and
%   maximum. For example,
%
%       F(X(1,1),X(1,2)) = Y(1)     and      F(X(2,1),X(2,2)) = Y(2)
%
% See also MAX2, MIN2, NORM.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = minandmax2@separableApprox(varargin{:});

end

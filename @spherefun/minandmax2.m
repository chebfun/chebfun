function varargout = minandmax2(varargin)
%MINANDMAX2   Find global minimum and maximum of a SPHEREFUN.
%   M = minandmax2(F) returns the minimum and maximum value of a SPHEREFUN
%   over its domain. M is a vector of length 2 such that 
%   M(1) = min(f(lambda,theta)) and M(2) = max(f(lambda,theta)).
%
%   [M, LOC] = minandmax2(F) also returns the position of the minimum and 
%   maximum. For example,
%
%       F(LOC(1,1),LOC(1,2)) = M(1)  and  F(LOC(2,1),LOC(2,2)) = M(2)
%
% See also SPHEREFUN/MAX2, SPHEREFUN/MIN2, SPHEREFUN/NORM.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = minandmax2@separableApprox(varargin{:});

end

function f = randnfun(dx, n, dom)
%RANDNFUN   Random smooth function
%   F = RANDNFUN(DX) returns a smooth chebfun on [-1,1] with maximum
%   wave number about 2pi/DX and standard normal distribution N(0,1)
%   at each point.  F can be regarded as one sample path of a Gaussian
%   process.  F is obtained by calling RANDNFUNTRIG on an interval
%   of about twice the length and restricting the result to [-1,1].
%
%   F = RANDNFUN(DX, N) returns a quasimatrix with N independent columns.
%
%   F = RANDNFUN(DX, N, DOM) returns a result with domain DOM = [A, B].
%
% Examples:
%
%   f = randnfun(0.1); std(f), plot(f)
%   plotcoeffs(f, '.'), xlim([0 200])
%
%   X = randnfun(.01,2); cov(X)
%
% See also RANDNFUNTRIG.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if nargin < 3
    dom = [-1 1];     % default domain
end

if nargin < 2
    n = 1;            % default number of columns
end

% Call RANDNFUNTRIG on interval of approximately double length.
% and then restrict the result to the prescribed interval.

m = 2*round(diff(dom)/dx)+1;
dom2 = dom(1) + [0 m*dx];
f = randnfuntrig(dx, n, dom2);
f = f{dom(1),dom(2)};

function f = randnfun2trig(dx, n, dom)
%RANDNFUN2TRIG   Random smooth periodic function
%   F = RANDNFUNTRIG(DX) returns a smooth periodic chebfun (trigfun)
%   on [-1,1] with maximum wave number about 2pi/DX and standard normal
%   distribution N(0,1) at each point.  F can be regarded as one sample
%   path of a Gaussian process.  It is obtained from a finite Fourier series
%   with independent normally distributed coefficients of equal variance.
%
%   F = RANDNFUNTRIG(DX, N) returns a quasimatrix with N independent columns.
%
%   F = RANDNFUNTRIG(DX, N, DOM) returns a result with domain DOM = [A, B].
%
% Examples:
%
%   f = randnfuntrig(0.1); std(f), plot(f)
%   plotcoeffs(f, '.'), xlim([-100 100])
%
%   X = randnfuntrig(.01,2); cov(X)
%
% See also RANDNFUN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if nargin < 3
    dom = [-1 1];     % default domain
end

if nargin < 2
    n = 1;            % default number of columns
end

% Although the output is real, complex arithmetic is used for the
% construction since the 'trig', 'coeffs' mode is only documented
% in this case.

m = round(diff(dom)/dx);
N = 2*m+1;
c = randn(N, N) + 1i*randn(N, N);
[x,y] = meshgrid(-m:m,-m:m);
c = c .* (x.^2+y.^2 <= m^2);
f = chebfun2(c/sqrt(pi*m^2/4), 'trig', 'coeffs');
f = real(f);

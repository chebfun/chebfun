function f = randnfuntrig(dx, n, dom)
%RANDNFUNTRIG   Random smooth periodic function
%   F = RANDNFUNTRIG(DX) returns a smooth periodic chebfun (trigfun)
%   on [-1,1] with maximum wave number about 2pi/DX and standard normal
%   distribution N(0,1) at each point.  F can be regarded as one sample
%   path of a Gaussian process.  It is obtained from a finite Fourier series
%   with independent normally distributed coefficients of equal variance.
%   (This code is preliminary and the definition may change in the future.)
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
c = randn(2*m+1, n) + 1i*randn(2*m+1, n);
c = (c + flipud(conj(c)))/2;
f = chebfun(c/sqrt(2*m+1), dom, 'trig', 'coeffs');

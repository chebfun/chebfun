function f = randnfuntrig(dt, dom)
%RANDNFUNTRIG   Random smooth periodic function
%   F = RANDNFUNTRIG(DT) returns a smooth periodic chebfun (trigfun)
%   on [-1,1] with space scale DT and standard normal distribution
%   N(0,1) at each point.  It can be regarded as a sample path of a
%   Gaussian process.  F is obtained from a finite Fourier series with
%   independent normally distributed coefficients of equal variance.
%
%   F = RANDNFUNTRIG(DT, DOM) returns a random periodic chebfun on
%   the interval DOM = [A, B].
%
% Example:
%   f = randnfuntrig(0.1);
%   std(f)
%   plot(f)
%   plotcoeffs(f, '.'), xlim([-100 100])
%
% See also RANDNFUN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if nargin == 1
    dom = [-1 1];
end

% Although the output is real, complex arithmetic is used for the
% construction since the 'trig', 'coeffs' mode is only documented
% in this case.

m = round(diff(dom)/dt);
c = randn(2*m+1, 1) + 1i*randn(2*m+1, 1);
c = (c + flipud(conj(c)))/2;
f = chebfun(c/sqrt(2*m+1), dom, 'trig', 'coeffs');

function f = randnfuntrig(varargin)
%RANDNFUNTRIG   Random smooth periodic function
%   F = RANDNFUNTRIG(DX) returns a smooth periodic chebfun (trigfun)
%   on [-1,1] with maximum wave number about 2pi/DX and standard normal
%   distribution N(0,1) at each point.  F can be regarded as one sample
%   path of a Gaussian process.  It is obtained from a finite Fourier series
%   with independent normally distributed coefficients of equal variance.
%   (This code is preliminary and the definition may change in the future.)
%
%   F = RANDNFUNTRIG(DX, DOM) returns a result with domain DOM = [A, B].
%
%   F = RANDNFUNTRIG(DX, N) returns a quasimatrix with N independent columns.
%
%   F = RANDNFUNTRIG(DX, 'norm') normalizes the output by dividing it
%   by SQRT(DX), so that white noise is approached in the limit DX -> 0.
%
%   F = RANDNFUNTRIG() uses the default value DX = 1.  Combinations such
%   as RANDNFUNTRIG(DOM), RANDNFUNTRIG('norm', DX) are allowed so long as
%   DX, if present, precedes N, if present.
%
% Examples:
%
%   f = randnfuntrig(0.1); std(f), plot(f)
%   plotcoeffs(f, '.'), xlim([-100 100])
%
%   X = randnfuntrig(.01,2); cov(X)
%
%   f = randnfuntrig(0.1,'norm',[0 10]); plot(cumsum(f))
%
% See also RANDNFUN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

[dx, n, dom, normalize] = parseInputs(varargin{:});

% Although the output is real, complex arithmetic is used for the
% construction since the 'trig', 'coeffs' mode is only documented
% in this case.

m = round(diff(dom)/dx);
c = randn(2*m+1, n) + 1i*randn(2*m+1, n);
c = (c + flipud(conj(c)))/2;
f = chebfun(c/sqrt(2*m+1), dom, 'trig', 'coeffs');
if normalize
    f = f/sqrt(dx);
end

end

function [dx, n, dom, normalize] = parseInputs(varargin)

dx = NaN;  
n = NaN;  
dom = NaN;
normalize = 0;

for j = 1:nargin
    v = varargin{j};
    if ischar(v)
        normalize = 1;
    elseif ~isscalar(v)
        dom = v;
    elseif isnan(dx) 
        dx = v;
    else
        n = v;
    end
end

if isnan(dx)
   dx = 1;          % default space scale
end
if isnan(n)
   n = 1;           % default number of columns
end
if isnan(dom)
   dom = [-1 1];    % default domain
end

end

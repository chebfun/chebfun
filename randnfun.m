function f = randnfun(varargin)
%RANDNFUN   Random smooth function
%   F = RANDNFUN(DX) returns a smooth chebfun on [-1,1] with maximum
%   wave number about 2pi/DX and standard normal distribution N(0,1)
%   at each point.  F can be regarded as one sample path of a Gaussian
%   process.  It is obtained by calling RANDNFUNTRIG on an interval
%   of length about 20% longer and restricting the result to [-1,1].
%
%   F = RANDNFUN(DX, DOM) returns a result with domain DOM = [A, B].
%
%   F = RANDNFUN(DX, N) returns a quasimatrix with N independent columns.
%
%   F = RANDNFUN(DX, 'norm') normalizes the output by dividing it
%   by SQRT(DX), so that white noise is approached in the limit DX -> 0.
%
%   F = RANDNFUN() uses the default value DX = 1.  Combinations such
%   as RANDNFUN(DOM), RANDNFUN('norm', DX) are allowed so long as
%   DX, if present, precedes N, if present.
%
% Examples:
%
%   f = randnfun(0.1); std(f), plot(f)
%   plotcoeffs(f, '.'), xlim([0 200])
%
%   X = randnfun(.01,2); cov(X)
%
%   f = randnfun(0.1,'norm',[0 10]); plot(cumsum(f))
%
% See also RANDNFUNTRIG.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

[dx, n, dom, normalize] = parseInputs(varargin{:});

% Call RANDNFUNTRIG on interval of approximately 20% greater length

m = 1.2*round(diff(dom)/dx)+1;
dom2 = dom(1) + [0 m*dx];
if normalize
    f = randnfuntrig(dx, n, dom2, 'norm');
else
    f = randnfuntrig(dx, n, dom2);
end

% Restrict the result to the prescribed interval.  Explicit rather than
% adaptive constructive gives speedup by a factor of about 3.

x = chebpts(5*m, dom);    % 5*m is large enough...
f = chebfun(f(x), dom);   % ...so that this is equivalent to f{dom(1),dom(2)}
f = simplify(f, 1e-13);   % loosened tolerance gives cleanest Chebyshev series

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


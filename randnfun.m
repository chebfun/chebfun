function f = randnfun(varargin)
%RANDNFUN   Random smooth function
%   F = RANDNFUN(DX) returns a smooth chebfun on [-1,1] with maximum
%   wave number about 2pi/DX and standard normal distribution N(0,1)
%   at each point.  F can be regarded as one sample path of a Gaussian
%   process.  It is obtained by calling RANDNFUNTRIG on an interval
%   of about twice the length and restricting the result to [-1,1].
%
%   F = RANDNFUN(DX, DOM) returns a result with domain DOM = [A, B].
%
%   F = RANDNFUN(DX, N) returns a quasimatrix with N independent columns.
%
%   Combinations RANDNFUN(DX, N, DOM) or RANDNFUN(DX, DOM, N) are
%   also allowed.  Commands RANDNFUN() or RANDNFUN(DOM) use the
%   default value DX = 1.
%
% Examples:
%
%   f = randnfun(0.1); std(f), plot(f)
%   plotcoeffs(f, '.'), xlim([0 200])
%
%   X = randnfun(.01,2); cov(X)
%
%   f = randnfun([0 100]); plot(cumsum(f))
%
% See also RANDNFUNTRIG.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

[dx, n, dom] = parseInputs(varargin{:});

% Call RANDNFUNTRIG on interval of approximately double length.
% and then restrict the result to the prescribed interval.

m = 2*round(diff(dom)/dx)+1;
dom2 = dom(1) + [0 m*dx];
f = randnfuntrig(dx, n, dom2);
f = f{dom(1),dom(2)};

end

function [dx, n, dom] = parseInputs(varargin)

dx = NaN;
n = NaN;
dom = NaN;

for j = 1:nargin
    v = varargin{j};
    if ~isscalar(v)
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


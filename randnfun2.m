function f = randnfun2(varargin)
%RANDNFUN2   Random smooth 2D function 
%   F = RANDNFUN2(DT) returns a smooth CHEBFUN2 on [-1,1,-1,1] with
%   maximum wave number about 2pi/DT in both the X and Y directions
%   and standard normal distribution N(0,1) at each point.  F is obtained
%   by calling RANDNFUN2TRIG on a domain of dimensions about 20% greater
%   and restricting the result to [-1,1,-1,1].    
%
%   F = RANDNFUN2(DT, DOM) returns a result with domain DOM = [A,B,C,D].
%
%   F = RANDNFUN2(DT, 'norm') normalizes the output by dividing it
%   by SQRT(DT), so that white noise is approached in the limit DT -> 0.
%
%   F = RANDNFUN2() uses the default value DT = 1.  Combinations such
%   as RANDNFUN2(DOM), RANDNFUN2('norm', DT) are allowed.
%
% Examples:
%
%   f = randnfun2(0.1); mean2(f), plot(f)
%   f = randnfun2(0.25); plot(roots(f)), axis equal
%
% See also RANDNFUN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

[dt, dom, normalize] = parseInputs(varargin{:});

% Call RANDNFUN2TRIG on domain of approximately 20% greater dimensions

if normalize
    f = randnfuntrig(dt2, n, dom, 'norm');
else
    f = randnfuntrig(dt2, n, dom);
end

% Restrict the result to the prescribed domain.  Explicit rather than
% adaptive construction probably gives speedup, but not tested yet.

x = chebpts2(
x = chebpts(5*m, dom);    % 5*m is large enough...
f = chebfun(f(x), dom);   % ...so that this is equivalent to f{dom(1),dom(2)}
f = simplify(f, 1e-13);   % loosened tolerance gives cleanest Chebyshev series



end

function [dt, dom, normalize] = parseInputs(varargin)

dt = NaN;  
dom = NaN;
normalize = 0;

for j = 1:nargin
    v = varargin{j};
    if ischar(v)
        normalize = 1;
    elseif ~isscalar(v)
        dom = v;
    else
        dt = v;
    end
end

if isnan(dt)
   dt = 1;               % default space scale
end
if isnan(dom)
   dom = [-1 1 -1 1];    % default domain
end

end

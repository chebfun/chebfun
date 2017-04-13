function f = randnfun2trig(varargin)
%RANDNFUN2TRIG   Random smooth periodic 2D function 
%   F = RANDNFUN2TRIG(DX) returns a smooth periodic chebfun2
%   on [-1,1,-1,1] with maximum wave number about 2pi/DX and standard normal
%   distribution N(0,1) at each point.  F is obtained from a finite 2D
%   Fourier series with independent normally distributed coefficients
%   of equal variance.
%
%   F = RANDNFUN2TRIG(DX, DOM) returns a result with domain DOM = [A,B,C,D].
%   **This is not yet implemented.**
%
%   F = RANDNFUN2TRIG(DX, 'norm') normalizes the output by dividing it
%   by SQRT(DX), so that white noise is approached in the limit DX -> 0.
%
%   F = RANDNFUN2TRIG() uses the default value DX = 1.  Combinations such
%   as RANDNFUN2TRIG(DOM), RANDNFUN2TRIG('norm', DX) are allowed.
%
% Examples:
%
%   f = randnfun2trig(0.1); mean2(f), plot(f)
%   f = randnfun2trig(0.25); plot(roots(f)), axis equal
%
% See also RANDNFUNTRIG.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

[dx, dom, normalize] = parseInputs(varargin{:});
if ~isequal(dom, [-1 1 -1 1])
    disp('for now, dom must be the unit square')
    dom = [-1 1 -1 1];
end

m = round(diff(dom(1:2))/dx);
mm = 2*m+1;
c = randn(mm, mm) + 1i*randn(mm, mm);  % random coefficients on a square
[x,y] = meshgrid(-m:m,-m:m);
c = c.*(x.^2 + y.^2 <= m^2);           % confine to a disk for isotropy
c = c/sqrt(nnz(c));                    % ensure var = 1 at each point
v = trigtech.coeffs2vals( ...          % 2D periodic coeffs to vals
    trigtech.coeffs2vals(c).').';
f = chebfun2(real(v), dom, 'trig');    % take real part, construct chebfun2
if normalize
    f = f/sqrt(dx);                    % normalize for 2D white noise
end

end

function [dx, dom, normalize] = parseInputs(varargin)

dx = NaN;  
dom = NaN;
normalize = 0;

for j = 1:nargin
    v = varargin{j};
    if ischar(v)
        normalize = 1;
    elseif ~isscalar(v)
        dom = v;
    else
        dx = v;
    end
end

if isnan(dx)
   dx = 1;               % default space scale
end
if isnan(dom)
   dom = [-1 1 -1 1];    % default domain
end

end

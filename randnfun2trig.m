function f = randnfun2trig(varargin)
%RANDNFUN2TRIG   Random smooth periodic 2D function 
%   F = RANDNFUN2TRIG(DT) returns a smooth periodic CHEBFUN2
%   on [-1,1,-1,1] with maximum wave number about 2pi/DT in both the
%   X and Y directions and standard normal distribution N(0,1) at
%   each point.  F is obtained from a finite 2D Fourier series with
%   independent normally distributed coefficients of equal variance.
%
%   F = RANDNFUN2TRIG(DT, DOM) returns a result with domain DOM = [A,B,C,D].
%
%   F = RANDNFUN2TRIG(DT, 'norm') normalizes the output by dividing it
%   by SQRT(DT), so that white noise is approached in the limit DT -> 0.
%
%   F = RANDNFUN2TRIG() uses the default value DT = 1.  Combinations such
%   as RANDNFUN2TRIG(DOM), RANDNFUN2TRIG('norm', DT) are allowed.
%
% Examples:
%
%   f = randnfun2trig(0.2); mean2(f), plot(f)
%   view(0,90), colormap(gray(2)), axis equal
%   plot(roots(f)), axis equal
%
% See also RANDNFUNTRIG.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

[dt, dom, normalize] = parseInputs(varargin{:});

m = round(diff(dom(1:2))/dt);
m2 = 2*m+1;
n = round(diff(dom(3:4))/dt);
n2 = 2*n+1;
c = randn(n2, m2) + 1i*randn(n2, m2);   % random coefficients on a square
[x,y] = meshgrid(-m:m,-n:n);
if ( m>0 & n>0 )
    c = c.*((x/m).^2 + (y/n).^2 <= 1);  % confine to a disk for isotropy
end
c = c/sqrt(nnz(c));                     % ensure var = 1 at each point
f = chebfun2(c, dom, 'coeffs', 'trig');
f = real(f);
if normalize
    f = f/sqrt(dt);                     % normalize for 2D white noise
end

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

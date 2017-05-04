function f = randnfun2(varargin)
%RANDNFUN2   Random smooth 2D function 
%   F = RANDNFUN2(DT) returns a smooth CHEBFUN2 on [-1,1,-1,1] with
%   maximum wave number about 2pi/DT in both the X and Y directions
%   and standard normal distribution N(0,1) at each point.  F is obtained
%   by calling RANDNFUN2(T, 'trig') on a domain of dimensions about 20%
%   greater and restricting the result to [-1,1,-1,1].    
%
%   F = RANDNFUN2(DT, DOM) returns a result with domain DOM = [A,B,C,D].
%
%   F = RANDNFUN2(DT, 'norm') normalizes the output by dividing it
%   by SQRT(DT).  (This may change.)
%
%   RANDNFUN(DT, 'trig') returns a random doubly periodic function.
%   This is defined by a finite bivariate Fourier series with independent
%   normally distributed coefficients of equal variance.
%
%   F = RANDNFUN2() uses the default value DT = 1.  Combinations such
%   as RANDNFUN2(DOM), RANDNFUN2('norm', DT) are allowed.
%
% Examples:
%
%   f = randnfun2(0.1); mean2(f), plot(f)
%   f = randnfun2(0.25); plot(roots(f)), axis equal
%
% See also RANDNFUN, RANDNFUNSPHERE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

[dt, dom, normalize, trig] = parseInputs(varargin{:});

if trig    % periodic case: random bivariate Fourier series

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

else       % nonperiodic case: call periodic case and restrict

    if dt == inf
        f = chebfun2(randn, dom);    % random constant function
        if normalize
            f = f/sqrt(dt);          % zero function
        end
        return
    end
    
    m = round(1.2*(dom(2)-dom(1))/dt+2);
    n = round(1.2*(dom(4)-dom(3))/dt+2);
    dom2 = [ dom(1) dom(1)+m*dt dom(3) dom(3)+n*dt ];

    if normalize
        f = randnfun2(dt, dom2, 'norm', 'trig');
    else
        f = randnfun2(dt, dom2, 'trig');
    end
    
           % restrict the result to the prescribed domain

    f = f{dom(1),dom(2),dom(3),dom(4)};

end
end

function [dt, dom, normalize, trig] = parseInputs(varargin)

dt = NaN;  
dom = NaN;
normalize = 0;
trig = 0;

for j = 1:nargin
    v = varargin{j};
    if ischar(v)
	if ( v(1) == 'n' )
            normalize = 1;
        elseif ( v(1) == 't' )
            trig = 1;
        else
            error('CHEBFUN:randnfun2','Unrecognized string input')
        end
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

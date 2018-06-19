function f = randnfun2(varargin)
%RANDNFUN2   Smooth random function in 2D
%   F = RANDNFUN2(LAMBDA) returns a CHEBFUN2 on [-1,1,-1,1] with
%   maximum frequency about 2pi/LAMBDA in both the X and Y directions
%   and standard normal distribution N(0,1) at each point.  F is obtained
%   by calling RANDNFUN2(T, 'trig') on a domain of dimensions about 20%
%   greater and restricting the result to [-1,1,-1,1].    
%
%   F = RANDNFUN2(LAMBDA, DOM) returns a result with domain DOM = [A,B,C,D].
%
%   F = RANDNFUN2(LAMBDA, 'big') normalizes the output by dividing it
%   by SQRT(LAMBDA).  (This may change.)
%
%   RANDNFUN(LAMBDA, 'trig') returns a random doubly periodic function.
%   This is defined by a finite bivariate Fourier series with independent
%   normally distributed coefficients of equal variance.
%
%   F = RANDNFUN2() uses the default value LAMBDA = 1.  Combinations such
%   as RANDNFUN2(DOM), RANDNFUN2('big', LAMBDA) are allowed.
%
% Examples:
%
%   f = randnfun2(0.1); mean2(f), plot(f)
%   f = randnfun2(0.25); plot(roots(f)), axis equal
%
% See also RANDNFUN, RANDNFUNSPHERE, RANDNFUNDISK.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

[lambda, dom, makebig, trig] = parseInputs(varargin{:});

if trig    % periodic case: random bivariate Fourier series

    m = round(diff(dom(1:2))/lambda);
    m2 = 2*m+1;
    n = round(diff(dom(3:4))/lambda);
    n2 = 2*n+1;
    c = randn(n2, m2) + 1i*randn(n2, m2);   % random coefficients on a square
    [x,y] = meshgrid(-m:m,-n:n);
    if ( m>0 & n>0 )
        c = c.*((x/m).^2 + (y/n).^2 <= 1);  % confine to a disk for isotropy
    end
    c = c/sqrt(nnz(c));                     % ensure var = 1 at each point
    f = chebfun2(c, dom, 'coeffs', 'trig');
    f = real(f);
    if makebig
        f = f/sqrt(lambda);                     % normalize for 2D white noise
    end

else       % nonperiodic case: call periodic case and restrict

    if lambda == inf
        f = chebfun2(randn, dom);    % random constant function
        if makebig
            f = f/sqrt(lambda);          % zero function
        end
        return
    end
    
    m = round(1.2*(dom(2)-dom(1))/lambda+2);
    n = round(1.2*(dom(4)-dom(3))/lambda+2);
    dom2 = [ dom(1) dom(1)+m*lambda dom(3) dom(3)+n*lambda ];

    if makebig
        f = randnfun2(lambda, dom2, 'big', 'trig');
    else
        f = randnfun2(lambda, dom2, 'trig');
    end
    
           % restrict the result to the prescribed domain

    f = f{dom(1),dom(2),dom(3),dom(4)};

end
end

function [lambda, dom, makebig, trig] = parseInputs(varargin)

lambda = NaN;  
dom = NaN;
makebig = 0;
trig = 0;

for j = 1:nargin
    v = varargin{j};
    if ischar(v)
        if ( v(1) == 'n' | v(1) == 'b')
            makebig = 1;
        elseif ( v(1) == 't' )
            trig = 1;
        else
            error('CHEBFUN:randnfun2','Unrecognized string input')
        end
    elseif ~isscalar(v)
        dom = v;
    else
        lambda = v;
    end
end

if isnan(lambda)
   lambda = 1;               % default space scale
end
if isnan(dom)
   dom = [-1 1 -1 1];    % default domain
end

end

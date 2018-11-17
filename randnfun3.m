function f = randnfun3(varargin)
%RANDNFUN3   Smooth random function in 3D
%   F = RANDNFUN3(LAMBDA) returns a CHEBFUN3 on [-1,1,-1,1,-1,1] with
%   maximum frequency about 2pi/LAMBDA in both the X, Y and Z directions
%   and standard normal distribution N(0,1) at each point.  F is obtained
%   by calling RANDNFUN3(T, 'trig') on a domain of dimensions about 20%
%   greater and restricting the result to [-1,1,-1,1,-1,1].    
%
% See also RANDNFUN, RANDNFUN2.

% Copyright 2018 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.


[lambda, dom, trig] = parseInputs(varargin{:});

if trig    % periodic case: random bivariate Fourier series

    m = round(diff(dom(1:2))/lambda);
    m2 = 2*m+1;
    n = round(diff(dom(3:4))/lambda);
    n2 = 2*n+1;
    p = round(diff(dom(5:6))/lambda);
    p2 = 2*p+1;
    c = randn(n2, m2, p2) + 1i*randn(n2, m2, p2);   % random coefficients on a cube
    [x,y,z] = ndgrid(-m:m,-n:n,-p:p);
    if ( m>0 && n>0 && p>0 )
        c = c.*((x/m).^2 + (y/n).^2 + (z/p).^2 <= 1);  % confine to a ball for isotropy
    end
    c = c/sqrt(nnz(c));                     % ensure var = 1 at each point
    f = chebfun3(c, dom, 'coeffs', 'trig');
    f = real(f);
    
else       % nonperiodic case: call periodic case and restrict

    if lambda == inf
        f = chebfun3(randn, dom);    % random constant function
        return
    end
    
    m = round(1.2*(dom(2)-dom(1))/lambda+2);
    n = round(1.2*(dom(4)-dom(3))/lambda+2);
    p = round(1.2*(dom(6)-dom(5))/lambda+2);
    dom2 = [ dom(1) dom(1)+m*lambda dom(3) dom(3)+n*lambda dom(5) dom(5)+p*lambda ];

    f = randnfun3(lambda, dom2, 'trig');
    
    % restrict the result to the prescribed domain

    f = f{dom(1),dom(2),dom(3),dom(4),dom(5),dom(6)};

end
end

function [lambda, dom, trig] = parseInputs(varargin)

lambda = NaN;  
dom = NaN;
trig = 0;

for j = 1:nargin
    v = varargin{j};
    if ischar(v)
        if ( v(1) == 't' )
            trig = 1;
        else
            error('CHEBFUN:randnfun3','Unrecognized string input')
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
   dom = [-1 1 -1 1 -1 1];    % default domain
end
end
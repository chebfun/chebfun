function r = diffrow(n, p, x, dom, type)
%DIFFROW     One row of spectral differentiation matrix.
%   DIFFROW(N, P, X) returns one row of DIFFMAT(N, P):
%         the first row if X = -1, the last row if X = +1.
%   DIFFROW(N, P, X, DOM) returns one row of DIFFMAT(N, P, DOM):
%         the first row if X = DOM(1), the last row if X = DOM(2).
%   DIFFROW(N, P, X, DOM, 'TRIG') returns the first row of
%         DIFFMAT(N, P, DOM, 'TRIG') if X = DOM(1).
%
% More general functionality, notably other values of X,
% has not yet been implemented.
%
% See also DIFFMAT, INTMAT, INTROW.

% This code has been written for Aurentz and Trefethen, 
% "Block operators and spectral discretizations", and is not
% yet fully compliant with standard Chebfun coding practices.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.


if ( nargin == 3 )                      % Chebyshev grid
    dom = [-1 1];
end

if ( nargin <= 4)                       % Chebyshev grid
    D = diffmat([2 n], p, dom);
    if ( x == dom(1) )
        r = D(1,:);
    elseif ( x == dom(2) )
        r = D(end,:);
    else
        error('CHEBFUN:diffow','illegal value of x')
    end
end

if ( nargin == 5)                       % trigonometric grid
    D = diffmat([n n], p, dom, 'trig');
    if ( x == dom(1) )
        r = D(1,:);
    else
        error('CHEBFUN:diffow','illegal value of x')
    end
end

end


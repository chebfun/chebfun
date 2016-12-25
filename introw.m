function r = introw(n, dom)
%INTROW     One row of spectral integration matrix.
%   INTROW(N) returns the last row of INTMAT(N), i.e., 
%         a row vector of N Clenshaw-Curtis quadrature coefficients.
%   INTROW(N, DOM) returns the last row of INTMAT(N, DOM), i.e.,
%         Clenshaw-Curtis coefficients scaled to the interval DOM.
%
% More general functionality, such as higher-order integrals
% or trigonometric quadrature, has not yet been implemented.
%
% See also DIFFMAT, DIFFROW, INTMAT.

% This code has been written for Aurentz and Trefethen, 
% "Block operators and spectral discretizations", and is not
% yet fully compliant with standard Chebfun coding practices.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 ) 
    dom = [-1 1];
end

[~,r] = chebpts(n,dom);

end


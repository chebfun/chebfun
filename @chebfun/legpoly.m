function out = legpoly(f, n)
%LEGPOLY    Legendre polynomial coefficients of a CHEBFUN.
%   A = LEGPOLY(F) returns the coefficients such that F_1 = A(1) P_N(x) + ... +
%   A(N) P_1(x) + A(N+1) P_0(x) where P_N(x) denotes the N-th Legendre
%   polynomial and F_1 denotes the first FUN of CHEBFUN F.
%
%   A = LEGPOLY(F, I) returns the coefficients for the I-th FUN.
%
%   There is also a LEGPOLY command in the Chebfun trunk directory, which
%   computes the CHEBFUN corresponding to the Legendre polynomial P_n(x).
%
% See also CHEBPOLY.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: Compute Legendre coefficients of a piecewise CHEBFUN by computing the
% required inner products.

% [TODO]: Make the second input be the number of returned coefficients, rather
% than the FUN of interest.

% Select a FUN:
if ( nargin == 1 )
    if ( numel(f.funs) > 1 )
        warning('CHEBFUN:legpoly:onefun', ...
               ['CHEBFUN has more than one FUN. ', ...
                'Returning Legendre coefficients for the first FUN only. ' ...
                'Use LEGPOLY(F, 1) to suppress this warning.'])
    end
    n = 1;
end

% Call FUN/LEGPOLY():
out = legpoly(f.funs{n});

end

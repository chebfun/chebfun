function out = legpoly(f, n)
%LEGPOLY    Legendre polynomial coefficients.
%   A = LEGPOLY(F) returns the coefficients such that F_1 = A(1) P_N(x) + ... +
%   A(N) P_1(x) + A(N+1) P_0(x) where P_N(x) denotes the N-th Legendre
%   polynomial and F_1 denotes the first fun of chebfun F.
%
%   A = LEGPOLY(F, I) returns the coefficients for the I-th fun.
%
%   There is also a LEGPOLY command in the Chebfun trunk directory, which
%   computes the chebfun corresponding to the Legendre polynomial P_n.
%
% See also CHEBPOLY.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Select a fun:
if ( nargin == 1 )
    if ( f.nfuns > 1 )
        warning('CHEBFUN:legpoly:onefun', ['Chebfun has more than one fun.', ...
                'Only the Legendre coefficients of the first one are returned.' ...
                 ' Use LEGPOLY(F, 1) to suppress this warning.'])
    end
    n = 1;
end

% Call FUN/LEGOPLY():
out = legpoly(f.funs{n});

end
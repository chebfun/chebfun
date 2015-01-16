function f = chebpoly(n, d, kind)
%CHEBPOLY   Chebyshev polynomial.
%   F = CHEBPOLY(N) returns the CHEBFUN corresponding to the Chebyshev
%   polynomials T_N(x) on [-1,1], where N may be a vector of nonnegative
%   integers.
%
%   F = CHEBPOLY(N, D), where D is an interval or a domain, gives the same
%   result scaled accordingly.
%
%   F = CHEBPOLY(N, KIND) or F = CHEBPOLY(N, D, KIND) returns Chebyshev
%   polynomials of the 1st kind, T_N(x)), when KIND = 1, and Chebyshev
%   polynomials of the 2nd kind, U_N(x)), when KIND = 2.
%
% See also LEGPOLY, FOURPOLY, and CHEBPTS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information. 

% TODO: This method requires a test.

% Defaults:
defaultKind = 1;

% Parse input
if ( any(n < 0) || any(mod(n, 1) ~= 0) )
    error('CHEBFUN:chebpoly:integern', ...
    'The first argument must be a vector of nonnegative integers.');
end

if ( nargin == 1 )
    d = chebfunpref().domain;
    kind = defaultKind;
elseif ( nargin == 2 )
    if ( numel(d) > 1 )
        kind = defaultKind;
    else
        kind = d;
        d = chebfunpref().domain;
    end
end    

% Cannot handle unbounded domains:
if ( any(isinf(d)) )
    error('CHEBFUN:chebpoly:infdomain', ...
    'Chebyshev polynomials are not defined over an unbounded domain.');
end

% [TODO]: 2nd-kind polynomials:
if ( kind == 2 )
    error('CHEBFUN:chebpoly:kind', ...
    '2nd-kind Chebyshev polynomials have not been implemented yet.');
end

% Construct the Chebyshev coefficients:
N = max(n)+1;
c = zeros(N, numel(n));
for k = 1:numel(n)
    c(n(k)+1, k) = 1;
end

% Construct a CHEBFUN from the coefficients:
f = chebfun(c, d([1, end]), 'coeffs');

% Introduce interior breakpoints:
if ( numel(d) > 2 )
    f = restrict(f, d);
end

% Transpose if required:
if ( size(n, 1) > 1 )
    f = f.'; 
end

end


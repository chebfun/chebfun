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
% See also LEGPOLY, TRIGPOLY, and CHEBPTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Add support for 3rd and 4th kind polynomials?

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

if ( ~( kind == 1 || kind == 2 ) )
    error('CHEBFUN:chebpoly:kind', ...
        'CHEBPOLY(N, KIND) only supports KIND = 1 or KIND = 2.');
end

% Cannot handle unbounded domains:
if ( any(isinf(d)) )
    error('CHEBFUN:chebpoly:infdomain', ...
    'Chebyshev polynomials are not defined over an unbounded domain.');
end

% Construct the Chebyshev coefficients:
N = max(n) + 1;
c = eye(N);
c = c(:,n+1);

% 2nd-kind polynomials:
if ( kind == 2 )
    % Use the recurrence relation from ultraS to map from Chebyshev U to T.
    c = ultraS.convertmat(N, 0, 0)\c;
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


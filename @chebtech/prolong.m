function f = prolong(f, nOut)
%PROLONG   Manually adjust the number of points used in a CHEBTECH.
%   G = PROLONG(F, N) returns a CHEBTECH G where LENGTH(G) = N and G represents
%   the same function as F but using more interpolation points/Chebyshev
%   coefficients then were stored in F.
%
%   If N < LENGTH(F) than the representation is compressed (by aliasing the
%   coefficients), which may result in loss of accuracy.
%
% See also CHEBTECH1/ALIAS, CHEBTECH2/ALIAS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Store the number of values the input function has:
nIn = size(f.coeffs, 1);

% nDiff is the number of new values needed (negative if compressing).
nDiff = nOut - nIn;

% Trivial case
if ( nDiff == 0 )
    % Nothing to do here!
    return
end

% Constant case:
if ( nIn == 1 )
    m = size(f.coeffs, 2);
    f.coeffs = [ f.coeffs(1,:) ; zeros(nDiff, m) ];
    return
end

% Prolong the coefficients using an FFT:
f.coeffs = f.alias(f.coeffs, nOut);

end
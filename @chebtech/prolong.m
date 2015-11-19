function f = prolong(f, nOut)
%PROLONG   Manually adjust the number of points used in a CHEBTECH.
%   G = PROLONG(F, N) returns a CHEBTECH G where LENGTH(G) = N and G represents
%   the same function as F but using more or less coefficients than F.
%
%   If N < LENGTH(F) the representation is compressed by chopping
%   coefficients, which may result in a loss of accuracy.
%
%   If N > LENGTH(F) the coefficients are padded with zeros.
%
% See also CHEBTECH1/ALIAS, CHEBTECH2/ALIAS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
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

% nDiff > 0
if ( nDiff > 0 )
    m = size(f.coeffs, 2);
    f.coeffs = [ f.coeffs ; zeros(nDiff, m) ];
    return
end

% nDiff < 0
if ( nDiff < 0 )
    m = max(nOut,0);
    f.coeffs = f.coeffs(1:m,:);;
    return
end

end

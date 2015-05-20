function out = iszero(f)
%ISZERO   True for zero DELTAFUN objects.
%   ISZERO(F) returns logical TRUE if F has a non-zero smooth part of some
%   non trivial delta functions

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Empty DELTAFUN is considered as a zero DELTAFUN:
if ( isempty(f) )
    out = 1;
    return
end

if ( ~iszero(f.funPart) )
    out = 0;
    return
end

if ( isempty(f.deltaLoc ) || isempty(f.deltaMag) )
    out = 1;
    return
end

pref = chebfunpref;
deltaTol = pref.deltaPrefs.deltaTol;
if ( max(abs(f.deltaMag(:))) < deltaTol )
    out = 1;
else
    out = 0;
end

end

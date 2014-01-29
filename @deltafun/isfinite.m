function out = isfinite(f)
%ISFINITE   Test if a DELTAFUN is bounded.
%   ISFINITE(F) returns FALSE if F has any non trivial delta functions or 
%   if the smooth part is infinite. It is TRUE otherwise.
%
% See also ISINF, ISNAN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check the smooth part:
if ( ~isempty(f.funPart) )
    if ( ~isfinite(f.funPart) )
        out = 0;
        return
    end
end

% Smooth part is finite, check the distributional part:
if ( isempty(f.deltaLoc) || isempty(f.deltaMag) )
    out = 1;
    return
end

% Get the tolerance:
pref = chebpref();
deltaTol = pref.deltaPrefs.deltaTol;

if ( max(abs(f.deltaMag)) < deltaTol )
    out = 1;
else
    out = 0;
end

end
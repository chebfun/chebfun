function out = anyDelta(f)
%ANYDELTA   True if a DELTAFUN object F has at least one delta function.
%   ANYDELTA(F) returns TRUE if the DELTAFUN object F has non-trivial delta
%   functions. The test is FALSE otherwise. This function uses the tolerance
%   provided by CHEBFUNPREF and uses that tolerance to decide whether a delta
%   function is trivial or not.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    out = 0;
    return 
end

% Get tolerance for deltas:
pref = chebfunpref();
deltaTol = pref.deltaPrefs.deltaTol;

% Check if f has no or only below tolerance delta functions.
out = any(abs(f.deltaMag(:)) > deltaTol);

end

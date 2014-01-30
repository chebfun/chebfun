function out = anyDelta(f)
%ANYDELTA   True if a DELTAFUN object F has atleast one delta function.
%   ANYDELTA(F) returns TRUE if the DELTAFUN object F has non-trivial delta
%   functions. The test is FALSE otherwise.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% TODO: NH: Comment on tolerance.

if ( isempty(f) )
    out = 0;
    return 
end

% Get tolerance for deltas:
pref = chebpref();
deltaTol = pref.deltaPrefs.deltaTol;

% Check if f has no or only below tolerance delta functions.
out = any(abs(f.deltaMag(:)) > deltaTol);

end
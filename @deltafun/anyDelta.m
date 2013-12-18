function out = anyDelta(f)
%ANYDELTA   True if a DELTAFUN object F has atleast one delta function.
%   ANYDELTA(F) returns TRUE if the DELTAFUN object F has non-trivial delta
%   functions. The test is FALSE otherwise.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Get tolerance for exponents:
tol = deltafun.pref.deltafun.deltaTol;

% A function is smooth if it has no or below tolerance delta functions.
out = any(abs(f.impulses(:)) > tol);
end
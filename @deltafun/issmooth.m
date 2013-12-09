function out = issmooth(f)
%ISEQUAL   True if a DELTAFUN object F has no delta functions.
%   ISSMOOTH(F) returns TRUE if the DELTAFUN object F has "negligible" delta
%   functions. The test is FALSE otherwise.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Get tolerance for exponents:
tol = deltafun.pref.deltafun.deltaTol;

% A function is smooth if it has no or below tolerance delta functions.
out = all(abs(f.delta.magnitude) < tol);
end
function seedRNG(s)
%SEEDRNG   Seed the MATLAB random number generator.
%   SEEDRNG(S) seeds the MATLAB random number generator with the integer value
%   S.  This function is meant to provide an interface for seeding the RNG that
%   is independent of the currently installed version of MATLAB.  It is used
%   solely for seeding the RNG in tests which use random numbers to ensure
%   deterministic behavior between runs.
%
%   Users should NOT call this function from their own code, as it may be
%   removed without warning in a future release.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( verLessThan('matlab', '7.12') )
    % Before R2011a.
    rand('seed', s);
    randn('seed', s);
else
    % R2011a and later
    rng('default');
    rng(s);
end

end

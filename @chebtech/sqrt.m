function g = sqrt(f)
%SQRT   Square root of a CHEBTECH.
%   SQRT(F) returns the square root of a CHEBTECH F.
%
%   This function assumes that the curve traced out by F in the complex plane
%   both (1) does not come too close to zero and (2) does not cross over the
%   branch cut in SQRT along the negative real axis.  That is, F should not
%   vanish at any point of [-1, 1], and the imaginary part of F should not
%   vanish at any point of (-1, 1) where the real part of F is negative.  If any
%   of these assumptions are violated, garbage may be returned with no warning.
%
% See also POWER.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% If F is complex and the imaginary part of F vanishes exactly at the domain 
% boundary, tiny rounding errors in evaluating the imaginary part at the
% boundary points can cause it to appear to change sign there.  If this happens,
% and if the real part of F is negative there, the endpoint function values will
% cross over the branch cut in SQRT, creating a jump discontinuity that should 
% not be there.  Enabling extrapolation resolves this issue by avoiding these 
% function evaluations:
pref = f.techPref();
if ( ~isreal(f) )
    pref.extrapolate = 1;
end

% Simply call the compose function:
g = compose(f, @sqrt, [], pref);

end
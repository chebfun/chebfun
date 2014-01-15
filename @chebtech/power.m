function g = power(f, b)
%.^   CHEBTECH power.
%   F.^G returns a CHEBTECH F to the scalar power G, a scalar F to the CHEBTECH
%   power G, or a CHEBTECH F to the CHEBTECH power G. F and or G may be complex.
%
%   This function assumes that the curve traced out by F in the complex plane
%   both (1) does not come too close to zero and (2) does not cross over the
%   branch cut in POWER along the negative real axis.  That is, F should not
%   vanish at any point of [-1, 1], and the imaginary part of F should not
%   vanish at any point of (-1, 1) where the real part of F is negative.  If any
%   of these assumptions are violated, garbage may be returned with no warning.
%
%   H = POWER(F, G) is called for the syntax 'F .^ G'.
%
% See also SQRT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% If B is a constant, cast it to a CHEBTECH:
if ( isnumeric(b) )
    b = f.make(@(x) 0*x+b, b, 1);
end

% If F is complex and the imaginary part of F vanishes exactly at the domain 
% boundary, tiny rounding errors in evaluating the imaginary part at the
% boundary points can cause it to appear to change sign there.  If this happens,
% and if the real part of F is negative there, the endpoint function values will
% cross over the branch cut in POWER, creating a jump discontinuity that should 
% not be there.  Enabling extrapolation resolves this issue by avoiding these 
% function evaluations:
pref = f.techPref();
if ( ~isreal(f) )
    pref.extrapolate = 1;
end

% Simply call the compose function:
g = compose(f, @power, b, pref);

end
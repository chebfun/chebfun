function g = power(f, b)
%.^   TRIGTECH power.
%   F.^G returns a TRIGTECH F to the scalar power G, a scalar F to the TRIGTECH
%   power G, or a TRIGTECH F to the TRIGTECH power G. F and or G may be complex.
%
%   This function assumes that the curve traced out by F in the complex plane
%   both (1) does not come too close to zero and (2) does not cross over the
%   branch cut in POWER along the negative real axis.  That is, F should not
%   vanish at any point of [-1, 1], and the imaginary part of F should not
%   vanish at any point of (-1, 1) where the real part of F is negative. If any
%   of these assumptions are violated, garbage may be returned with no warning.
%
%   H = POWER(F, G) is called for the syntax 'F .^ G'.
%
% See also SQRT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

pref = f.techPref();

% Simply call the compose function:
if ( isa(f, 'trigtech') && isa(b, 'trigtech') )
    % Both F and G are TRIGTECHs:
	g = compose(f, @power, b, pref);
elseif ( isa(f, 'trigtech') )
    % F is TRIGTECH and G is constant: 
	g = compose(f, @(f) power(f, b), [], pref);
else
    % F is constant and G is TRIGTECH: 
	g = compose(b, @(b) power(f, b), [], pref);
end

end

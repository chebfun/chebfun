function p = angle(f, pref)
%ANGLE   Chebfun phase angle.
%   ANGLE(H) returns the phase angles, in radians, of a complex-valued CHEBFUN.
%
% See also ABS, UNWRAP, ATAN2.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 2 )
    pref = chebfunpref();
end

% Call ATAN2():
p = atan2(imag(f), real(f), pref);

end

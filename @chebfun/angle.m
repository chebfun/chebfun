function p = angle(f)
%ANGLE  Chebfun phase angle.
%   ANGLE(H) returns the phase angles, in radians, of a complex-valued chebfun.
%
% See also ANGLE, CHEBFUN/ABS, CHEBFUN/UNWRAP.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

p = atan2(imag(f), real(f));

end
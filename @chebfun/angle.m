function p = angle(f)
%ANGLE  Chebfun phase angle.
%   ANGLE(H) returns the phase angles, in radians, of a complex-valued CHEBFUN.
%
% See also ABS, UNWRAP, ATAN2.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% [TODO]: This requires ATAN2()!
p = atan2(imag(f), real(f));

end
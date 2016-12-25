function f = rotate(h, alpha)
%ROTATE   Rotate a DISKFUN by alpha, where alpha is an angle in radians.
%   Y = ROTATE(H, ALPHA) rotates H by the angle ALPHA. Positive ALPHA
%   rotates H in the counterclockwise direction and negative ALPHA rotates 
%   H in the clockwise direction.
% 
% See also DISKFUN/CIRCSHIFT.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call CIRCSHIFT: 
if ( nargin == 1 ) 
    alpha = 0; 
end

f = circshift(h, alpha);

end
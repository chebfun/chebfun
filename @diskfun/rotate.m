function f = rotate(h, alpha)
%ROTATE: rotates a DISKFUN by alpha, where alpha is an angle in radians.
%
%   Y = ROTATE(H, ALPHA) rotates H by the angle ALPHA. Positive ALPHA
%   rotates H in the counterclockwise direction, negative ALPHA rotates H
%   in the clockwise direction.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


if nargin == 1
    return
else   
    fr = @(t,r) feval(h, t-alpha, r, 'polar'); 
    f = diskfun(fr, 'polar');
end
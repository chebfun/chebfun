function f = circshift(h, alpha)
%CIRCSHIFT: Circular (periodic) shift of a diskfun by alpha units in the 
%   angular (polar) direction.
%
%   f = CIRCSHIFT(h, ALPHA) rotates h by the angle ALPHA. Positive ALPHA
%   rotates H in the counterclockwise direction, negative ALPHA rotates H
%   in the clockwise direction.
%
%   See also DISKFUN/ROTATE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


if nargin == 1
    f = h;
    return
else   
    f = h;  
    %get trigtech
    rtechs = f.rows.funs{1}.onefun;
    rtechs = circshift(rtechs, alpha/pi);
    f.rows.funs{1}.onefun = rtechs;
end
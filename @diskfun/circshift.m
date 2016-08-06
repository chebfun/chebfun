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

% This is the only way to get this to work for now:
f = diskfun(@(t, r) feval(h, t-alpha, r), 'polar'); %rotate by alpha

% It is faster to make it work like this, but this is breaking encapsulation: 
% if nargin == 1
%     f = h;
%     return
% else   
%     f = h;  
%     %get trigtech
%     rtechs = f.rows.funs{1}.onefun;
%     rtechs = circshift(rtechs, alpha/pi);
%     f.rows.funs{1}.onefun = rtechs;
%     r = f.rows; 
%     r.pointValues = chebfun.getValuesAtBreakpoints(r.funs, r.domain);
%     f.rows = r; 
% end

end
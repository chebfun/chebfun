function f = circshift(h, alpha)
%CIRCSHIFT   Circular (periodic) shift of a diskfun by alpha units in the 
%   angular (polar) direction.
%
%   f = CIRCSHIFT(h, ALPHA) rotates h by the angle ALPHA. Positive ALPHA
%   rotates H in the counterclockwise direction, negative ALPHA rotates H
%   in the clockwise direction.
%
%   See also DISKFUN/ROTATE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% This is what we want to do.  However, we can do it much master by
% drilling down to the trigtechs representing the rows.  If one is worried
% about encapsulation then use this commented out code segement and remove
% the code starting at the if.
% h = diskfun(@(t, r) feval(h, t-alpha, r, 'polar'), 'polar'); %rotate by alpha

if nargin == 1
    f = h;
    return
else   
    f = h;  
    % Operate at the trigtech level to make things fast
    rtechs = circshift(f.rows.funs{1}.onefun, alpha/pi);
    f.rows.funs{1}.onefun = rtechs;
    % Weird feval behavior in CHEBFUN requires this:
    f.rows.pointValues = feval(rtechs, [-1;1]); 
end

end
function f = fliplr(g)
%FLIPLR   Flip/reverse a DISKFUN over the x axis. 
%
%   G = FLIPLR( F ) returns a DISKFUN G that is flipped over the y-axis
%   that is, G(x,y) = F(-x,y).
%
% See also DISKFUN/FLIPUD, DISKFUN/ROTATE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%check for empty

if ( isempty( g ) ) 
    return
end 

f = diskfun(@(x,y) feval(g, -x, y, 'cart'));


end

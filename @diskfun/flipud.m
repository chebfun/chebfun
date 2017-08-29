function f = flipud(g)
%FLIPUD   Flip/reverse a DISKFUN over the x-axis.
%
%   G = FLIPUD(F) returns a DISKFUN G that is flipped over the 
%   x-axis, that is G(x,y) = F(x, -y).
%
% See also DISKFUN/FLIPLR, DISKFUN/FLIPDIM, DISKFUN/ROTATE. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty( g ) ) 
    f = diskfun;
    return
end 
f = g; 
f.rows = flipud(g.rows); 

end

function g = flipxy(f)
% FLIPXY   G = FLIPXY(F) returns the diskfun G(Y, X) given F(X,Y). 
% 
% See also FLIPUD, FLIPLR, and FLIPDIM

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
if ( isempty(f) )
    g = diskfun;
    return
end

g = rotate(f, -pi/4);
g.rows = flipud(g.rows); 
g = rotate(g, pi/4); 

end

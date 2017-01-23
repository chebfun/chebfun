function grid = getGrid(~, N, dom)
%GETGRID   Returns a grid correspoding to a SPINOP3 object.
%
% See also SPINOP3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

xx = trigpts(N, dom(1:2));
yy = trigpts(N, dom(3:4));
zz = trigpts(N, dom(5:6));
[xx, yy, zz] = meshgrid(xx, yy, zz);
grid = {xx; yy; zz};

end
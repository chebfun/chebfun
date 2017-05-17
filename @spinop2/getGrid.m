function grid = getGrid(~, N, dom)
%GETGRID   Returns a grid correspoding to a SPINOP2 object.
%
% See also SPINOP2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

xx = trigpts(N, dom(1:2));
yy = trigpts(N, dom(3:4));
[xx, yy] = meshgrid(xx, yy);
grid = {xx; yy};
    
end
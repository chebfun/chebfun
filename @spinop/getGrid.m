function grid = getGrid(~, N, dom)
%GETGRID   Returns a grid correspoding to a SPINOP object.
%
% See also SPINOP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

xx = trigpts(N, dom(1:2));
grid = {xx};
    
end
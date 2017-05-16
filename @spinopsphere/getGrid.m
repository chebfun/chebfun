function grid = getGrid(~, N, dom)
%GETGRID   Returns a grid correspoding to a SPINOPSPHERE object.
%
% See also SPINOPSPHERE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note: SPINOPSPHERE and SPHEREFUN objects live on [-pi,pi]x[0,pi] so the domain
% DOM that is passed is DOM = [-pi,pi]x[0,pi]. However, the grid we use for the 
% computation is a doubled-up grid corresponding to [-pi,pi]^2.

ll = trigpts(N, dom(1:2)); % [-pi, pi]
tt = trigpts(N, [-pi pi]); % impose [-pi pi] and not dom(3:4) = [0 pi]
[ll, tt] = meshgrid(ll, tt);
grid = {ll; tt};
    
end
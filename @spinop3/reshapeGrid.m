function grid = reshapeGrid(~, grid)
%RESHAPEGRID   Add the repeated endpoints to a 2D periodic grid.

% Note: We use periodic discretizations in 1D/2D/3D. These discretizatioins 
% do not include the repeated endpoints. Before plotting data with 
% SPINOPERATOR/INITIALIZEMOVIE, we add these repeated endpoints to the grid. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the points:
xx = grid{1};
yy = grid{2};
zz = grid{3};

% Add the endpoints:
xx = [xx, 2*xx(:,end,:) - xx(:,end-1,:)];
xx =  [xx; xx(1,:,:)];
xx = cat(3, xx, xx(:,:,1));
yy = [yy; 2*yy(end,:,:) - yy(end-1,:,:)];
yy = [yy, yy(:,1,:)];
yy = cat(3, yy, yy(:,:,1));
zz = cat(3, zz, 2*zz(:,:,end) - zz(:,:,end-1));
zz = [zz; zz(1,:,:)];
zz = [zz, zz(:,1,:)];

% Output the new grid:
grid{1} = xx;
grid{2} = yy;
grid{3} = zz;

end
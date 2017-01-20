function grid = reshapeGrid(~, grid)
%RESHAPEGRID   Add the repeated endpoints to a 2D periodic grid.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note: We use periodic discretizations in 1D/2D/3D. These discretizatioins 
% do not include the repeated endpoints. Before plotting data with 
% SPINOPERATOR/INITIALIZEMOVIE, we add these repeated endpoints to the grid. 

% Get the points:
xx = grid{1};
yy = grid{2};

% Add the endpoints:
xx = [xx, 2*xx(:,end) - xx(:,end-1)];
xx = [xx; xx(1,:)];
yy = [yy; 2*yy(end,:) - yy(end-1,:)];
yy = [yy, yy(:,1)];

% Output the new grid:
grid{1} = xx;
grid{2} = yy;

end
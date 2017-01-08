function grid = reshapeGrid(~, grid)
%RESHAPEGRID   Add the repeated endpoint to a 1D periodic grid.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note: We use periodic discretizations in 1D/2D/3D. These discretizatioins 
% do not include the repeated endpoints. Before plotting data with 
% SPINOPERATOR/INITIALIZEMOVIE, we add these repeated endpoints to the grid. 

% Get the points:
xx = grid{1};

% Add the endpoint:
xx = [xx; 2*xx(end) - xx(end-1)];

% Output the new grid:
grid{1} = xx;

end
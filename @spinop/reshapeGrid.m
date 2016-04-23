function grid = reshapeGrid(~, grid)
%RESHAPEGRID   Add the repeated endpoint to a 1D periodic grid.

% Get the points:
xx = grid{1};

% Add the endpoint:
xx = [xx; 2*xx(end) - xx(end-1)];

% Output the new grid:
grid{1} = xx;

end
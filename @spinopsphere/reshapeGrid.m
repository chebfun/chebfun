function grid = reshapeGrid(~, grid)
%RESHAPEGRID    Add the repeated endpoints in the lambda- and theta-directions 
%and extract half of the points in the theta-direction.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note: We use periodic discretizations in longitude (lambda) and latitude
% (theta). These dicretizations do not include the repeated endpoints. Before 
% plotting data with SPINOPERATOR/INITIALIZEMOVIE, we add these repeated 
% endpoints to the grid, and extract half of the points in theta, those that
% correspond to [0, pi]. (We use the DFS method with doubled-up [-pi, pi] 
% theta-grid.)

% Get the points:
ll = grid{1};
tt = grid{2};

% Get the number of grid points:
N = length(ll);

% Add the endpoints in lambda and theta:
ll = [ll, 2*ll(:,end) - ll(:,end-1)];
ll = [ll; ll(1,:)];
tt = [tt; 2*tt(end,:) - tt(end-1,:)];
tt = [tt, tt(:,1)];

% Extract half of the points in the theta-direction:
ll = ll(N/2+1:end,:);
tt = tt(N/2+1:end,:);

% Output the new grid:
grid{1} = ll;
grid{2} = tt;

end
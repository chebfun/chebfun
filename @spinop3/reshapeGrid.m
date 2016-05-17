function grid = reshapeGrid(~, grid)
%RESHAPEGRID   Add the repeated endpoints to a 2D periodic grid.

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
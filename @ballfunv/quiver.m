function quiver(v,varargin)
% QUIVER Plot of a BALLFUNV on the ball

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Take a vector field Vx, Vy, Vz and plot it
[Vx,Vy,Vz] = v.comp{:};

[m,n,p] = size(v);

% Create the grid
r   = chebpts(m);
lam = pi*trigpts(n);
th  = [pi*trigpts(p); pi];

% Convert the coeffs of the vector field to vals
Vx = ballfun.coeffs2vals(Vx.coeffs);
Vy = ballfun.coeffs2vals(Vy.coeffs);
Vz = ballfun.coeffs2vals(Vz.coeffs);

% Remove the doubles-up data: [0,1]x[-pi,pi[x[0,pi]
Vx = real(Vx(floor(m/2)+1:end,:,[ceil(p/2)+1:end 1]));
Vy = real(Vy(floor(m/2)+1:end,:,[ceil(p/2)+1:end 1]));
Vz = real(Vz(floor(m/2)+1:end,:,[ceil(p/2)+1:end 1]));

% Resize the indices and grid
r = r(floor(m/2)+1:end);
th = th(ceil(p/2)+1:end);

% Convert the spherical points to cartesian points
[Lam,R,Th] = meshgrid(lam,r,th);
x = R.*cos(Lam).*sin(Th);
y = R.*sin(Lam).*sin(Th);
z = R.*cos(Th);

q = quiver3(x,y,z,Vx,Vy,Vz, 'AutoScaleFactor',0.5);

% Color the vectors according to their magnitude
if (nargin==1 || (nargin==2 && varargin{1}~="Color") ||  (nargin>2 && varargin{1}~="Color" && varargin{2}~="Color"))
    % Compute the magnitude of the vectors
    mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), ...
                reshape(q.WData, numel(q.UData), [])).^2, 2));
            
    % Scale the colormap to data
    caxis([0,max(mags)]);
    colorbar;
            
    % Get the current colormap
    currentColormap = colormap(jet);
     
    % Now determine the color to make each arrow using a colormap
    % The colors scale to the axis of the colorbar
    clims = num2cell(get(gca, 'clim'));
    [~, ~, ind] = histcounts(mags, linspace(clims{:}, size(currentColormap, 1)));

    % Now map this to a colormap to get RGB
    cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
    cmap(:,:,4) = 255;
    cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

    % We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
    set(q.Head, ...
        'ColorBinding', 'interpolated', ...
        'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'

    % We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
    set(q.Tail, ...
        'ColorBinding', 'interpolated', ...
        'ColorData', reshape(cmap(1:2,:,:), [], 4).');
end
hold on;

% Sphere grid
if (nargin==1 || (nargin==2 && varargin{1}~="Grid") ||  (nargin>2 && varargin{1}~="Grid" && varargin{2}~="Grid"))
    [x,y,z] = sphere(16);
    surf(x, y, z, 'linestyle',':','FaceAlpha', 0);
end

% Add label
xlabel('X')
ylabel('Y')
zlabel('Z')

% Axis
axis([-1 1 -1 1 -1 1]);
daspect([1 1 1]);
hold off;
end

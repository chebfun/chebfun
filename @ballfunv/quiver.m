function varargout = quiver(v, varargin)
%QUIVER   Quiver plot of BALLFUNV.
%   QUIVER(F) plots the vector velocity field of F. QUIVER automatically
%   attempts to scale the arrows to fit within the grid.
%
%   QUIVER(F,S) automatically scales the arrows to fit within the grid and then
%   stretches them by S.  Use S=0 to plot the arrows without the automatic
%   scaling. The arrows are on a uniform grid.
%
%   QUIVER(...,LINESPEC) uses the plot linestyle specified for the velocity
%   vectors.  Any marker in LINESPEC is drawn at the base instead of an arrow on
%   the tip.  Use a marker of '.' to specify no marker at all.  Use color 
%   specifiers to specify the color of the arrows.  See PLOT for other possibilities.
%
%   QUIVER(...,'numpts',N) plots arrows on a N by N by N uniform grid.
%
% See also BALLFUN/PLOT.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check if one component is empty
if isempty(v)
    isVempty = 1;
else
    isVempty = cellfun(@isempty, v.comp, 'UniformOutput', false);
    isVempty = any(cell2mat(isVempty));
end
if isVempty
    error('CHEBFUN:BALLFUNV:quiver:isempty','ballfunv must not have an empty component.');
end

% Default parameters
numpts = 25;
scale = 2.5;

% Number of points to plot
j = 1;
argin = {};
color = true;

while ( ~isempty( varargin ) )
    if strcmpi( varargin{1}, 'numpts' )
        numpts = varargin{2};
        varargin(1:2) = [];
    elseif strcmpi( varargin{1}, 'color' )
        color = false;
        argin{j} = varargin{1};
        argin{j+1} = varargin{2};
        varargin(1:2) = [];
        j = j+2;
    elseif ismember( varargin{1}, ['r','g','b','c','m','y','k','w'])
        color = false;
        argin{j} = varargin{1};
        varargin(1) = [];
        j = j+1;
    elseif isnumeric( varargin{1} )
        scale = varargin{1};
        varargin(1) = [];
    else
        argin{j} = varargin{1};
        varargin(1) = [];
        j = j+1;
    end
end
varargin = {scale, argin{:}};

% Create the grid
m = numpts; n = numpts; p = numpts;
r = chebpts(m);

% Resize the indices and grid
r = r(floor(m/2)+1:end);
m = length(r);

% Create the vector in cartesian system
Vxx = [];
Vyy = [];
Vzz = [];

% Initialize the array of equally spaced points
xx = [];
yy = [];
zz = [];

% Find the equally spaced points in the ball
for i = 1:m
    % Coordinates of the points on the sphere of radius r(i)
    Nth = max(ceil(p*r(i)/2),1);
    th_i  = linspace(0,pi,Nth);
    
    for k = 1:Nth
        Dth = min(th_i(k),abs(th_i(k)-pi));
        Nlam = max(ceil(n*r(i)*Dth*2/pi),1);
        lam_i = trigpts(Nlam)*pi;
        
        % Get the values of the vector field at these points
        [VX, VY, VZ] = fevalm(v,r(i),lam_i,th_i(k));
        
        Vxx = [Vxx;VX(1,:,1).'];
        Vyy = [Vyy;VY(1,:,1).'];
        Vzz = [Vzz;VZ(1,:,1).'];
        x = r(i)*cos(lam_i)*sin(th_i(k));
        y = r(i)*sin(lam_i).*sin(th_i(k));
        z = repmat(r(i)*cos(th_i(k)),Nlam,1);
        xx = [xx;x];
        yy = [yy;y];
        zz = [zz;z];
    end
end

q = quiver3(xx,yy,zz,real(Vxx),real(Vyy),real(Vzz), varargin{:});

% Color the vectors according to their magnitude
if  color
    % Compute the magnitude of the vectors
    mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), ...
                reshape(q.WData, numel(q.UData), [])).^2, 2));
            
    % Scale the colormap to data
    caxis([0,max(mags)]);
            
    % Get the current colormap
    currentColormap = colormap();
     
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

% % Axis
axis([-1 1 -1 1 -1 1]);
daspect([1 1 1]);
axis square

if ( nargout > 0 )
    varargout = { q }; 
end

end
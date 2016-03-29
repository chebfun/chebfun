function varargout = plotEarth(varargin)
%PLOTEARTH  Plots the landmasses of earth.
%   PLOTEARTH plots outlines of the landmasses of earth on a sphere using
%   black solid lines.
%
%   PLOTEARTH(LINESPEC) uses the uses the color and marker from the 
%   line specification string 'LineSpec' (See PLOT for possibilities)
%
% See also spherefun/surf, spherefun/plot, spherefun/contour

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Land masses are stored in the data file CoastData.mat
try
    x = load('CoastData.mat', 'coast');
catch
    error('CHEBFUN:SPHEREFUN:PLOTEARTH:coastDataNotFound',...
        ['File containing the coast line data could not be found. '...
        'Try reinstalling chebfun.']);
end

if ( nargin > 0 )
    linespec = varargin{:};
else
    linespec = 'k-';
end

% Get the hold state of the current axis:
holdState = ishold;

if ( ~holdState )
    % Generate a unit sphere.
    [XX,YY,ZZ] = sphere(101);
    
    % Color of the sphere will be white:
    clr = [255 255 255]/255;
    
    % Plot the sphere, make it slightly smaller than unit so lines
    % show up more clearly.
    scl = 0.99;
    
    surf(scl*XX, scl*YY, scl*ZZ, 1+0*XX, 'EdgeColor', 'None', ...
        'FaceColor', clr);
    daspect([1 1 1]);
    hold on
end

h = plot3(x.coast(:,1), x.coast(:,2), x.coast(:,3), linespec);

if ( ~holdState )
    hold off
end

if ( nargout > 0 )
    varargout = { h }; 
end

end
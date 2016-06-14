function varargout = plot(f, varargin)
%PLOT   Plots a CHEBFUN3 object.
%   Plot(F) creates six contour plots at boundary cross sections 
%   of the domain of F.
% 
%   If F is complex-valued, then six phase portraits are plotted at 
%   boundary cross sections.
%
%   Example: 
%   plot(chebfun3(@(x,y,z) sin(i*x+z)+cos(y))); campos([-10 -11 -8])
%
% See also CHEBFUN3/SLICE, CHEBFUN3/SCAN, CHEBFUN3/ISOSURFACE, and 
% CHEBFUN3/SURF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

holdState = ishold;
dom = f.domain;
numpts = 151;
[xx, yy, zz] = meshgrid(linspace(dom(1), dom(2), numpts), ...
    linspace(dom(3), dom(4), numpts), linspace(dom(5), dom(6), numpts));
v = feval(f, xx, yy, zz);

% Slices are at the bounds of the cube:
xSlices = [dom(1) dom(2)];
ySlices = [dom(3) dom(4)];
zSlices = [dom(5) dom(6)];

if ( isreal(v) )
    h = slice(xx, yy, zz, v, xSlices, ySlices, zSlices); 
    shading interp
    colorbar
    axis(gca, 'tight')
else
    % Phase portraits if F is complex-valued:
    h = slice(xx, yy, zz, angle(-v), xSlices, ySlices, zSlices);
    set(h, 'EdgeColor','none')
    caxis([-pi pi])
    colormap('hsv')
    axis('equal') 
end

if ( ~holdState )
    hold off
end
if ( nargout > 0 )
    varargout = {h};
end

end
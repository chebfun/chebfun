function varargout = plot(f, varargin)
%surf3(x,y,z,varargin)
%PLOT   Plots a CHEBFUN3 object.
%   Plot(F) creates six contour plots at boundary cross sections 
%   of the domain of F.
% 
%   If F is complex-valued, then plot(F) plots six phase portraits at 
%   boundary cross sections.
%
%   PLOT(X,Y,Z,F) where X, Y, Z, and F are all CHEBFUN3 objects, plots F as
%   above but also allows for domains that are more general than the cube.
%
%   PLOT(X,Y,Z) puts F = Z, i.e., colors the plot according to the
%   Z-values. This is analogous to MATLAB's surf(x,y,z).
%
%   Example: 
%   plot(chebfun3(@(x,y,z) sin(i*x+z)+cos(y))); campos([-10 -11 -8])
%
% See also CHEBFUN3/SLICE, CHEBFUN3/SCAN, CHEBFUN3/ISOSURFACE, and 
% CHEBFUN3/SURF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

holdState = ishold;

if nargin <= 2
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
    return
end

if nargin == 3
    x = f;
    y = varargin{1};
    z = varargin{2};
elseif nargin == 4
    x = f;
    y = varargin{1};
    z = varargin{2};
    f = varargin{3};
elseif nargin >= 5
    x = f;
    y = varargin{1};
    z = varargin{2};
    f = varargin{3};
end

d = x.domain;
index.type = '()';
index.subs = {[d(1)]  ':'  ':'};
surf(subsref(x,index),subsref(y,index),subsref(z,index),subsref(f,index));
hold on
index.subs = {[d(2)]  ':'  ':'};
surf(subsref(x,index),subsref(y,index),subsref(z,index),subsref(f,index));

index.subs = {':' [d(3)]  ':'};
surf(subsref(x,index),subsref(y,index),subsref(z,index),subsref(f,index));
index.subs = {':' [d(4)]  ':'};
surf(subsref(x,index),subsref(y,index),subsref(z,index),subsref(f,index));

index.subs = {':'  ':' [d(5)]};
surf(subsref(x,index),subsref(y,index),subsref(z,index),subsref(f,index));
index.subs = {':'  ':' [d(6)]};
h = surf(subsref(x,index),subsref(y,index),subsref(z,index),subsref(f,index));
hold off
    
if ( ~holdState )
    hold off
end
if ( nargout > 0 )
    varargout = {h};
end

end
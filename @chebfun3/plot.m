function varargout = plot(f, varargin)
%PLOT   Makes six contour plots forming boundary of the domain.
% 
%   If f is complex-valued, then six phase portraits are plotted.
%
%   Example: 
%   plot(chebfun3(@(x,y,z) sin(i*x+z)+cos(y))); campos([-10 -11 -8])
%
%   See also CHEBFUN3/SLICE, SCAN, CHEBFUN3/ISOSURFACE, and CHEBFUN3/SURF.

holdState = ishold;
dom = f.domain;
numpts = 151;
[xx, yy, zz] = meshgrid(linspace(dom(1), dom(2), numpts), ...
    linspace(dom(3), dom(4),numpts), linspace(dom(5), dom(6), numpts));
v = feval(f, xx, yy, zz);
% if ( ~isreal(v) )
%     dicideWhatToDo
%     v = abs(v); % complex magnitude
% end

xSlices = [dom(1) dom(2)];
ySlices = [dom(3) dom(4)];
zSlices = [dom(5) dom(6)];
if ( isreal(v) )
    h = slice(xx, yy, zz, v, xSlices, ySlices, zSlices); 
    shading interp
    colorbar
    % set(gca, 'xticklabel',{[]}) 
    % set(gca, 'yticklabel',{[]}) 
    % set(gca, 'zticklabel',{[]}) 
    axis(gca, 'tight')
else
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
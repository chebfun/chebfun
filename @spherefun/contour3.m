function varargout = contour3( f, varargin )
%CONTOUR3   3-D contour plot of a SPHEREFUN.
%   CONTOUR3(F) is a contour plot of F treating the values of F as heights
%   above the sphere, analagous to SURF(F, 'projection', 'bumpy'). A
%   contour plot shows the level curves of F for some values V. The values
%   V are chosen automatically.
%
%   CONTOUR3(F, N) draws N contour lines, overriding the automatic number.
%   The values V are still chosen automatically.
%   
%   CONTOUR3(F, V) draws LENGTH(V) contour lines at the values specified in
%   the vector V. Use CONTOUR3(F, [V V]) to compute a single contour at the
%   level V.
%
%   CONTOUR3(F, 'NUMPTS', N) plots the contour lines on an N by N uniform
%   grid. If NUMPTS is not given then we plot on a 200 by 200 grid.
%
% See also CONTOUR, CONTOURF.

% Copyright 2020 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )  % Empty check.
    contour3([]);
    return
end

holdState = ishold;

% Minimum number of plotting points:
minplotnum = 200;

% Extract from the inputs the user defined options: 
j = 1; 
argin = {};
while ( ~isempty(varargin) )
    if ( strcmpi(varargin{1}, 'numpts') ) % If given numpts then use them.
        minplotnum = varargin{2};
        varargin(1:2) = [];
    else
        argin{j} = varargin{1};
        varargin(1) = [];
        j = j+1;
    end
end

if ( isa(f, 'spherefun') )
    dom = f.domain;
    % Evaluate at equally spaced grid: 
    x = linspace(dom(1), dom(2), minplotnum);
    y = linspace(dom(3), dom(4), minplotnum);
    vals = sample(f, minplotnum-1, minplotnum);
    vals = [vals vals(:,1)];
else
    error('CHEBFUN:SPHEREFUN:contour3:inputs', ...
        'Input must be a spherefun.');
end

if ( iscolat(f) )
    xc = @(ll,tt,rr) rr.*cos(ll).*sin(tt);
    yc = @(ll,tt,rr) rr.*sin(ll).*sin(tt);
    zc = @(tt,rr) rr.*cos(tt);
else
    xc = @(ll,tt,rr) rr.*cos(ll).*cos(tt);
    yc = @(ll,tt,rr) rr.*sin(ll).*cos(tt);
    zc = @(tt,rr) rr.*sin(tt);
end

% Use contour rather than contourc so that it can handle parsing the inputs
% correctly.
[c, h] = contour( x', y', vals, argin{:} );

% Extract out the options we need to plot the contours with plot3.
LW = 'LineWidth'; 
lw = h.LineWidth;
LS = 'LineStyle'; 
ls = h.LineStyle;
LC = 'Color'; 
lc = h.LineColor;
levelList = h.LevelList;
clrmap = parula(numel(levelList));

% Remove the contour plot that was generated.
delete(h);

scl = 0.15;           % Match the bumpy parameter in SPHEREFUN/SURF().
lim = [-1-scl 1+scl]; % Pad the axis limits.
m = minandmax2(f);

% If the plot is not being added to another, then set the axis properties.
if ( ~holdState )
    xlim(lim), ylim(lim), zlim(lim)
    daspect([1 1 1])
    grid on, box off
    hold on
end
view(3)

cl = size(c,2);
k = 1;
while ( k < cl )
    kl = c(2,k);
    v = k+1:k+kl;

    % Bump out the radial components according to the function values.
    rr = rescale(c(1,k), 1-scl, 1+scl, 'InputMin', m(1), 'InputMax', m(2));
    xv = xc(c(1,v), c(2,v), rr);
    yv = yc(c(1,v), c(2,v), rr);
    zv = zc(c(2,v), rr);

    % If the line color is a float then we are plotting all contours in a
    % single color.
    if ( isfloat(lc) )
        plot3(xv, yv, zv, LW, lw, LC, lc, LS, ls);
    else
        % We need to plot each contour in a color using the default 
        % colormap.
        % Determine the color for the level being plotted.
        clr = clrmap(abs(c(1, k) - levelList) < 10*eps, :);
        plot3(xv, yv, zv, LW, lw, LC, clr, LS, ls);
    end

    k = k+kl+1;
end

if ( ~holdState )
    hold off
end

% Return plot handle if appropriate.
if ( nargout >= 1 )
    warning('CHEBFUN:SPHEREFUN:contour3:outputs', ...
        'Outputs from contour3 are not supported');
    varargout = { [], [] };
end

end

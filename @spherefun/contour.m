function varargout = contour(f, varargin)
%CONTOUR    contour plot of a SPHEREFUN.
%   CONTOUR(F) is a contour plot of F treating the values of F as heights
%   above or below the unit sphere. Contours are the level curves of F for 
%   some values V. The values V are chosen automatically.
%
%   CONTOUR(F, N) draws N contour lines, overriding the automatic number. 
%   The values V are still chosen automatically.
%   
%   CONTOUR(F, V) draws LENGTH(V) contour lines at the values specified in 
%   vector V. Use contour(F, [v, v]) to compute a single contour at the 
%   level v.
%
%   CONTOUR(F, 'NUMPTS', N) plots the contour lines using samples of F on
%   an N-by-N grid in intrinsic coordinates. If NUMPTS is not given then we
%   plot on a 200-by-200 grid.
%
%   CONTOUR(F, 'PIVOTS', STR) plots the contour lines with the pivot
%   locations used during constructor. STR should indicate the marker to be
%   used for the pivots (see PLOT for options).
%
% See also PLOT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )  % Empty check.
    contour([]);
    return
end

holdState = ishold;

% Minimum number of plotting points:
minplotnum = 200;
doPivotPlot = 0;

% Extract from the inputs the user defined options: 
j = 1; 
argin = {};
while ( ~isempty(varargin) )
    if ( strcmpi(varargin{1}, 'numpts') ) % If given numpts then use them.
        minplotnum = varargin{2};
        varargin(1:2) = [];
    elseif ( strcmpi(varargin{1}, 'pivots') ) % If given numpts then use them.
        doPivotPlot = 1;
        if ( length(varargin) < 2 ) 
            error('CHEBFUN:SPHEREFUN:contour:pivotStyle', ...
                'Pivot style undefined.')
        end
        argin{j} = varargin{2};
        varargin(1:2) = [];
    else
        argin{j} = varargin{1};
        varargin(1) = [];
        j = j+1;
    end
end

% Did the user want a plot of the pivot locations?
if ( doPivotPlot )    % Do pivot plot. 
    if ( (~isempty(argin)) && (length(argin{1}) < 5) )
        
        % Column, row, pivot plot
        plot(f, argin{:}), 
        hold on
        argin(1) = [];
        contour(f, argin{:})
        return
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
    error('CHEBFUN:SPHEREFUN:contour:inputs', ...
        'Input must be a spherefun.');
    
end

if ( iscolat(f) )
    xc = @(ll,tt) cos(ll).*sin(tt);
    yc = @(ll,tt) sin(ll).*sin(tt);
    zc = @(tt) cos(tt);
else
    xc = @(ll,tt) cos(ll).*cos(tt);
    yc = @(ll,tt) sin(ll).*cos(tt);
    zc = @(tt) sin(tt);
end

% Use contour rather than contourc so that it can handle parsing the inputs
% correctly.
[C,H] = contour( x', y', vals, argin{:} );

% Extract out the options we need to plot the contours with plot3.
LW = 'LineWidth'; 
lw = H.LineWidth;
LS = 'LineStyle'; 
ls = H.LineStyle;
LC = 'Color'; 
lc = H.LineColor;
levelList = H.LevelList;
clrmap = parula(numel(levelList));

% Remove the contour plot that was generated.
delete(H);

% If the plot is not being added to another then plot a solid 
% sphere so the lines are more easily discernable.
if ( ~holdState )
    
    % Generate a unit sphere.
    [XX,YY,ZZ] = sphere(101);
    
    % Color of the sphere will be off-white:
    clr = [250 250 250]/255;
    
    % Plot the sphere, make it slightly smaller than unit so lines
    % show up more clearly.
    scl = 0.99;
    surf(scl*XX, scl*YY, scl*ZZ, 1+0*XX, 'EdgeColor', 'None', ...
        'FaceColor', clr);
    daspect([1 1 1]);
    hold on
end

cl = size(C,2);
k = 1;
while ( k < cl )
    kl = C(2,k);
    v = k+1:k+kl;
    xv = xc(C(1, v), C(2, v));
    yv = yc(C(1, v), C(2, v));
    zv = zc(C(2, v));
    
    % If the line color is a float then we are plotting all contours in a
    % single color.
    if ( isfloat(lc) )
        plot3(xv,yv,zv,LW,lw,LC,lc,LS,ls);
    else
        % We need to plot each contour in a color using the default 
        % colormap.
        % Determine the color for the level being plotted.
        clr = clrmap(abs(C(1, k) - levelList) < 10*eps, :);
        plot3(xv, yv, zv, LW, lw, LC, clr, LS, ls);
    end        
    k = k+kl+1;
end
if ( ~holdState )
    hold off
end

% return plot handle if appropriate.
if ( nargout >= 1 )
    warning('CHEBFUN:SPHEREFUN:contour:outputs', ...
        'Outputs from contour are not supported');
    varargout = { [], [] };
end

end
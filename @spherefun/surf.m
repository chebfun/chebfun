function varargout = surf(f, varargin)
%SURF  Surface plot of a SPHEREFUN.
%   SURF(F) plots F on the surface of a sphere.
%
%   SURF(..., 'GRID', S) plots F on the surface of a sphere and includes a
%   grid with standard line properties specifed by the string S.
%
%   SURF(..., 'PROJECTION', type) plots F using one of the following map
%   projections specified by the string type:
%     'sphere'          : 3D plot on the sphere (default)
%     'bumpy'           : 3D plot on the sphere with bumps
%     'equirectangular' : Standard latitude-longitude cylindrical projection
%     'hammer'          : Hammer (or Hammer-Aitoff) pseudoazimuthal projection
%     'albers'          : Albers conic projection
%     'eckert2'         : Eckert II pseudocylindrical projection
%     'winkel3'         : Winkel III (or tripel) pseudoazimuthal projection
%     'sinusoidal'      : Sinusoidal (or Mercator equal-area) pseudocylindrical 
%                         projection
%
%   SURF(..., 'NUMPOINTS', N) plots F using an N-by-N latitude-longitude
%   sample grid.
%
%   SURF(..., 'PropertyName', PropertyValue,...) sets the value of the
%   specified surface property. Multiple property values can be set with a
%   single statement. See surf for more details.
%
%   H = SURF(...) returns a handle to a surface plot object.
%
%   Note: to plot columns and rows of a SPHEREFUN, use PLOT(..., '.-').
%
%   Examples:
%      f = cheb.gallerysphere('vortices');
%      surf(f, 'grid')
%
%      f = spherefun.sphharm(4,0) + sqrt(5/7)*spherefun.sphharm(4,4);
%      surf(f, 'projection', 'hammer', 'grid')
%
%      f = spherefun.sphharm(10,4);
%      surf(f,'projection','bumpy'), camlight, axis off
%
% See also MESH, PLOT.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Empty check:
if ( isempty(f) )
    h = surf([]);
    if ( nargout == 1 )
        varargout = { h };
    end
    return
end

% Is the plot currently being held?
plotOnHold = ishold;
% How dense to make the samples?
minPlotNum = 200;
% How dense to make the grid lines?
numGridPtsLam = 24;
numGridPtsTh = 12;
% Default grid lines:
gridLineType = 'k-';
% Default is to not add a grid
addGrid = 0;
% Default plotting is on the sphere
plotType = 'sphere';
% Default plotting options
defaultOpts = {'facecolor', 'interp', 'edgealpha', .5, 'edgecolor', 'none'};

% Number of points to plot
j = 1;
argin = {};
while ( ~isempty(varargin) )
    if strcmpi(varargin{1}, 'numpts')
        minPlotNum = varargin{2};
        varargin(1:2) = [];
    elseif strcmpi(varargin{1}, 'grid')
        addGrid = 1;
        if mod(nargin, 2) == 0 
        % e.g., surf(f, 'grid') or surf(f, 'grid', 'projection', 'hammer')                       
            gridLineType = 'k-';
            varargin(1) = [];
        else
        % e.g. surf(f, 'grid', 'k-')
            gridLineType = varargin{2};
            varargin(1:2) = [];
        end
    elseif strcmpi(varargin{1}, 'projection')
        plotType = varargin{2};
        if ( ~any(strcmpi(plotType, {'sphere', 'bumpy'})) )
            % Only add grid automatically for 2D projections.
            addGrid = 1;
        end
        varargin(1:2) = [];
    else
        argin{j} = varargin{1};
        varargin(1) = [];
        j = j+1;
    end
end

if ( isempty(argin) )
    argin = {};
end

if ( isa(f,'spherefun') )
    % Default now is everything uses colatitude.
    l = linspace(-pi, pi, minPlotNum);
    t = linspace(0, pi, minPlotNum);
    C = fevalm(f, l, t);
    
    % Longitude-colatitude evaluation points
    [ll, tt] = meshgrid(l, t);
    
    % Make some corrections to C for prettier plotting.
    if ( norm(C - C(1,1),inf) < 1e-10 )
        % If vals are very close up to round off then the color scale is
        % hugely distorted. This fixes that.
        [n, m] = size(C);
        C = C(1,1)*ones(n, m);
    end
    
    % Meshgrid for plotting lines of longitude.
    [llgl, ttgl] = meshgrid(linspace(-pi, pi, numGridPtsLam+1),...
        linspace(0, pi, minPlotNum));
    % Meshgrid for plotting lines of latitude.
    [llgt, ttgt] = meshgrid(linspace(-pi, pi, minPlotNum),...
        linspace(0, pi, numGridPtsTh+1));
    
    
    % Plot on the surface of the sphere: 3D plot
    if any(strcmpi(plotType, {'sphere', 'bumpy'}))
        vv = ones(size(ll));
        lim = [-1 1];
        if strcmpi(plotType, 'bumpy')
            % Bump out the radial components according to function values.
            scl = 0.15;
            vv = vv + rescale(C, -scl, scl); % Only bump out by scl.
            lim = lim + [-scl scl];          % Pad the axis limits.
        end
        [xx,yy,zz] = sph2cart(ll,pi/2 - tt,vv);
        h = surf(xx, yy, zz, C, defaultOpts{:}, argin{:});
        if addGrid == 1
            [xxg,yyg,zzg] = sph2cart(llgl,pi/2-ttgl,1+0*llgl);
            hold on;
            plot3(xxg,yyg,zzg,gridLineType);
            [xxg,yyg,zzg] = sph2cart(llgt,pi/2-ttgt,1+0*llgt);
            plot3(xxg',yyg',zzg',gridLineType);
            if ~plotOnHold
                hold off;
            end
        end
        xlim(lim)
        ylim(lim)
        zlim(lim)
    else
        % Compute the mapped points
        [xh,yh] = sph2map(plotType,ll,tt);
        h = surf(xh, yh, 0*xh, C, defaultOpts{:}, argin{:});
        if addGrid == 1
            [xg,yg] = sph2map(plotType,llgl,ttgl);
            hold on;
            plot(xg,yg,gridLineType);
            [xg,yg] = sph2map(plotType,llgt',ttgt');
            plot(xg,yg,gridLineType);
            if ~plotOnHold
                hold off;
            end
        end
        if ~plotOnHold
            view(2)
            axis('tight','off');
        end
    end
    % Make the aspect ratio equal if the plot is not currently being
    % held.
    if ( ~plotOnHold )
        daspect([1 1 1]);
    end
    %     else
    %         % Pass this along to the surf function in separableApprox.
    %         h = surf@separableApprox(f, varargin{:});
    %     end
end

if ( nargout > 0 )
    varargout = { h };
end

end

function [xh,yh] = sph2map(type,lam,th)
%SPH2MAP 2D map projection from standar spherical coordinates.
%
% [XH,YH] = sph2map(TYPE,LAM,TH) maps points on the surface of the sphere
%   represented in spherical coordiantes to points in the 2D coordinate
%   system defined by the type.  Options for type are
%   'HAMMER'
%   'ALBERS'
%
%   Note that LAM is the longitude (in radians) of the points and TH is
%   colatitude (in radians), i.e. 0 <= TH <= pi, of the points.

type = lower(type);
switch type
    case 'equirectangular'
        xh = lam;
        yh = th;
    case 'hammer'
        xh = 2*sqrt(2)*sin(th).*sin(lam/2)./sqrt(1 + sin(th).*cos(lam/2));
        yh = sqrt(2)*cos(th)./sqrt(1 + sin(th).*cos(lam/2));
    case 'albers'
        th0 = pi/2;
        lam0 = 0;
        th1 = pi/2;
        th2 = pi/6;
        n = 0.5*(cos(th1) + cos(th2));
        phi = n*(lam - lam0);
        C = sin(th1)^2 + 2*n*cos(th1)^2;
        rho = sqrt(C - 2*n*cos(th))/n;
        rho0 = sqrt(C - 2*n*cos(th0))/n;
        xh = rho.*sin(phi);
        yh = rho0 - rho.*cos(phi);
    case 'eckert2'
        th = pi/2-th;
        xh = 2*lam.*sqrt((4 - 3*sin(abs(th)))/6*pi);
        yh = sign(th).*(sqrt(2*pi/3)*(2 - sqrt(4 - 3*sin(abs(th)))));
    case 'winkel3'
        th = pi/2-th;
        th1 = acos(2/pi);
        sincalpha = mysinc(acos(cos(th).*cos(lam/2)));
        xh = 0.5*(lam*cos(th1) + 2*cos(th).*sin(lam/2)./sincalpha);
        yh = 0.5*(th + sin(th)./sincalpha);
    case 'mercator'
        th = pi/2-th;
        xh = lam;
        yh = log(tan(pi/4 + th/2));
    case 'sinusoidal'
        th = pi/2-th;
        xh = lam.*cos(th);
        yh = th;
    otherwise
        xh = lam;
        yh = th;
end
end

function out = mysinc(x)
% Deal with the removable singularity at 0 explicitly.
out = sin(x)./x;
out(x == 0) = 1;
end

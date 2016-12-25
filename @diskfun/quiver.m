function varargout = quiver( F, G, varargin )
%QUIVER   Quiver plot of DISKFUN.
%   QUIVER(F,G) plots the vector velocity field of (F,G). QUIVER automatically
%   attempts to scale the arrows to fit within the grid. This returns the 
%   same plot as QUIVER([F ; G]).
%
%   QUIVER(F,G, S) automatically scales the arrows to fit within the grid and then
%   stretches them by S.  Use S=0 to plot the arrows without the automatic
%   scaling. 
%
%   QUIVER(X,Y,F,...) is the same as QUIVER(F,...) except the arrows are on the
%   grid given in X and Y.
%
%   QUIVER(...,LINESPEC) uses the plot linestyle specified for the velocity
%   vectors.  Any marker in LINESPEC is drawn at the base instead of an arrow on
%   the tip.  Use a marker of '.' to specify no marker at all.  See PLOT for
%   other possibilities.
%
%   QUIVER(...,'numpts',N) plots arrows on a N by N grid.
%
%   H = QUIVER(...) returns a quivergroup handle.
%
%   See also, DISKFUNV/QUIVER

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

numpts = 30;

% Empty check:
if ( isempty( F ) )
    quiver( [] )
    return
end

if ( isempty( varargin ) )
    varargin = {};
end

holdState = ishold;
if ( ~holdState )
    hold on;
end

if ( ~holdState )
    % Generate a unit disk
    N = 200;
    th = linspace(-pi, pi, N)';
    r = exp( 1i*th );
    plot(real(r), imag(r), 'k--');
end

% Number of points to plot
j = 1;
argin = {};
while ( ~isempty( varargin ) )
    if strcmpi( varargin{1}, 'numpts' )
        numpts = varargin{2};
        varargin(1:2) = [];
    else
        argin{j} = varargin{1};
        varargin(1) = [];
        j = j+1;
    end
end
varargin = argin;

if ( isa(F, 'diskfun') && isa(G, 'diskfun' ) )            % quiver(F,G,...)
    
        % Plot quiver with arrows at equally spaced points:
        [xx, yy] = diskpts(numpts);
        vals1 = feval(F, xx, yy, 'cart');
        vals2 = feval(G, xx, yy, 'cart');
        h = quiver(xx, yy, vals1,vals2, varargin{:});
        if ( ~holdState )
            axis tight;            
            axis( max([max(abs(xlim)) abs(ylim)])*[-1 1 -1 1] );
            axis equal;
        end
    
elseif ( nargin >= 3 )                 % quiver(x,y,F,...)
    
    % First two arguments contain arrow locations: we assume these are
    % Cartesian coords. 
    xx = F;
    yy = G;
    
    if ( isa(varargin{3}, 'diskfun') && isa(varargin{4},'diskfun')  )
        
        F = varargin{3};
        G = varargin{4};
            vals1 = feval(F, xx, yy, 'cart');
            vals2 = feval(G, xx, yy, 'cart');
            h = quiver( xx, yy, vals1, vals2, varargin{3:end} );
            if ( ~holdState )
                axis tight;
                axis(max([max(abs(xlim)) abs(ylim)])*[-1 1 -1 1]);
                axis equal;
            end
            
    else
        
        error('DISKFUN:DISKFUNV:quiver:inputs', ...
                                  'Third and fourth arguments should be diskfuns.');
        
    end
    
end

if ( ~holdState )
    hold off;
end

if ( nargout > 0 )
    varargout = {h};
end

end

% Generates a nice set of points on the unit disk that are roughly equally
% spaced.
function [xx,yy] = diskpts(numpts)
% Idea is to use a polar grid (r,theta), but instead of using the same
% number of points in theta for every r, we make it a function of r.  We
% start with the origin, then as r increases from the origin there will be
% 6, 9, 15, 21, ... until the outer r=1 radius is reached.  This appears to
% give a nice sampling of the disk.

% The number of radii to use to get approximately numpts is 
n = floor( numpts/sqrt(3) );

% Increment for the radii
dr = 1/n;

% Add the origin
xx = 0; 
yy = 0;

% Handle the second ring from the origin a bit differently (i.e. do 6
% points instead of 3).
th = trigpts(6,[-pi pi]);
xx = [xx; dr*cos(th)];
yy = [yy; dr*sin(th)];

% Add points 3*(2*k-1) points for radii k.
for k = 2:n
    th = trigpts(3*(2*k-1), [-pi pi]);
    xx = [xx ; dr*k*cos(th)];
    yy = [yy ; dr*k*sin(th)];
end
end

% ANOTHER OPTION, BUT I DON'T THINK IT'S AS NICE.
% % Generates a nice set of points on the unit disk using the technique
% % described in Section 3.3 of 
% % D. Calhoun, C. Helzel, R. J. LeVeque. Logically rectangular grids and
% % finite volume methods for PDEs in circular and spherical domains. SIAM
% % Review Vol. 50, Issue 4, pp. 723-752. (2008)
% function [xx,yy] = diskpts(numpts)
% 
% % Generate equally spaced points over [-1,1]x[-1,1].  These will be mapped
% % to the unit disk according the the algorithm given in Section 3.3. of the
% % paper referenced above.
% [xc, yc] = meshgrid(linspace(-1,1,numpts));
% 
% r1 = 1;   % map [-1,1] x [-1,1] to circle of radius r1
% d = max(abs(xc),abs(yc));
% r = sqrt(xc.^2 + yc.^2);
% r = max(r,1e-10);
% xx = r1 * d .* xc./r;
% yy = r1 * d .* yc./r;
% w = d.^2;
% xx = w.*xx + (1-w).*xc/sqrt(2);
% yy = w.*yy + (1-w).*yc/sqrt(2);
% 
% end
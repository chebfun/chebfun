function varargout = quiver( F, varargin )
%QUIVER   Quiver plot of DISKFUNV.
%   QUIVER(F) plots the vector velocity field of F. QUIVER automatically
%   attempts to scale the arrows to fit within the grid. The arrows are on a
%   uniform grid.
%
% BELOW HERE NOT YET IMPLEMENTED:
%   QUIVER(F,S) automatically scales the arrows to fit within the grid and then
%   stretches them by S.  Use S=0 to plot the arrows without the automatic
%   scaling. The arrows are on a uniform grid.
%
%   QUIVER(X,Y,F,...) is the same as QUIVER(F,...) except the arrows are on the
%   grid given in X and Y.
%
%   QUIVER(...,LINESPEC) uses the plot linestyle specified for the velocity
%   vectors.  Any marker in LINESPEC is drawn at the base instead of an arrow on
%   the tip.  Use a marker of '.' to specify no marker at all.  See PLOT for
%   other possibilities.
%
%   QUIVER(...,'numpts',N) plots arrows on a N by N uniform grid.
%
%   H = QUIVER(...) returns a quivergroup handle.
%
%   If F is a DISKFUN with three non-zero components then this calls
%   QUIVER3. (this is currently not implemented in diskfun)
%
% See also QUIVER3.
% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

numpts = 30;

% Empty check:
if ( isempty( F ) )
    quiver([])
    return
end

if ( isempty(varargin) )
    varargin = {};
end

holdState = ishold;
if ~holdState
    hold on;
end

if ~holdState
    %
    % Generate a unit disk
    N = 200;
    th = linspace(-pi,pi,N)';
    r = exp(1i*th);
    plot(real(r),imag(r), 'k--');
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

if ( isa(F, 'diskfunv') )             % quiver(F,...)
    
    nF = F.nComponents;
    if ( nF == 3 )
        h = quiver3(F, varargin{:});   % Call quiver3 instead.
    else
        % Plot quiver with arrows at equally spaced points:
        [xx, yy] = diskpts(numpts);
        F1 = F.components{1}; F2 = F.components{2};
        vals1 = feval(F1, xx, yy, 'cart');
        vals2 = feval(F2, xx, yy, 'cart');
        h = quiver(xx, yy, vals1,vals2, varargin{:});
        if ~holdState
            axis tight;            
            axis(max([max(abs(xlim)) abs(ylim)])*[-1 1 -1 1]);
            axis equal;
        end
    end
    
elseif ( nargin >= 3 )                 % quiver(x,y,F,...)
    
    % First two arguments contain arrow locations: %for diskfun we need to
    % know if this is cartesian or polar...for now we assume cartesian
    
    xx = F;
    yy = varargin{1};
    
    if ( isa(varargin{2}, 'diskfunv') )
        F = varargin{2};
        nF = F.nComponents;
        if ( nF == 3 )
            h = quiver3(F,varargin{:}); % Call quiver3 instead.
        else
            F1 = F.components{1}; F2 = F.components{2};
            vals1 = feval(F1, xx, yy, 'cart');
            vals2 = feval(F2, xx, yy, 'cart');
            h = quiver( xx, yy, vals1, vals2, varargin{3:end} );
            if ~holdState
                axis tight;
                axis(max([max(abs(xlim)) abs(ylim)])*[-1 1 -1 1]);
                axis equal;
            end
        end
    else
        error('DISKFUN:DISKFUNV:quiver:inputs', ...
            'Third argument should be a diskfunv.');
    end
    
end

if ~holdState
    hold off;
end

if ( nargout > 0 )
    varargout = {h};
end

end

% Generates a nice set of points on the unit disk using the technique
% described in Section 3.3 of 
% D. Calhoun, C. Helzel, R. J. LeVeque. Logically rectangular grids and
% finite volume methods for PDEs in circular and spherical domains. SIAM
% Review Vol. 50, Issue 4, pp. 723-752. (2008)
function [xx,yy] = diskpts(numpts)

% Generate equally spaced points over [-1,1]x[-1,1].  These will be mapped
% to the unit disk according the the algorithm given in Section 3.3. of the
% paper referenced above.
[xc, yc] = meshgrid(linspace(-1,1,numpts));

r1 = 1;   % map [-1,1] x [-1,1] to circle of radius r1
d = max(abs(xc),abs(yc));
r = sqrt(xc.^2 + yc.^2);
r = max(r,1e-10);
xx = r1 * d .* xc./r;
yy = r1 * d .* yc./r;
w = d.^2;
xx = w.*xx + (1-w).*xc/sqrt(2);
yy = w.*yy + (1-w).*yc/sqrt(2);

end

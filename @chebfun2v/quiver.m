function varargout = quiver( F, varargin )
%QUIVER   Quiver plot of CHEBFUN2V.
%   QUIVER(F) plots the vector velocity field of F. QUIVER automatically
%   attempts to scale the arrows to fit within the grid. The arrows are on a
%   uniform grid.
%
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
%   If F is a CHEBFUN2V with three non-zero components then this calls QUIVER3.
%
% See also QUIVER3.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

numpts = 10;

% Empty check:
if ( isempty( F ) )
    quiver([])
    return
end

if ( isempty(varargin) )
    varargin = {};
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

if ( isa(F, 'chebfun2v') )             % quiver(F,...)
    
    nF = F.nComponents;
    if ( nF == 3 )
        h = quiver3(F, varargin{:});   % Call quiver3 instead.
    else
        % Plot quiver with arrows at equally spaced points:
        dom = F.components{1}.domain;
        x = linspace(dom(1), dom(2), numpts);
        y = linspace(dom(3), dom(4), numpts);
        [xx, yy] = meshgrid(x, y);
        h = quiver(xx, yy, F, varargin{:});
    end
    
elseif ( nargin >= 3 )                 % quiver(x,y,F,...)
    
    % First two arguments contain arrow locations:
    xx = F;
    yy = varargin{1};
    
    if ( isa(varargin{2}, 'chebfun2v') )
        F = varargin{2};
        nF = F.nComponents;
        if ( nF == 3 )
            h = quiver3(F,varargin{:}); % Call quiver3 instead.
        else
            F1 = F.components{1}; F2 = F.components{2};
            vals1 = feval(F1, xx, yy);
            vals2 = feval(F2, xx, yy);
            dom = F1.domain;
            h = quiver( xx, yy, vals1, vals2, varargin{3:end} );
            axis(1.1*dom);
        end
    else
        error('CHEBFUN:CHEBFUN2V:quiver:inputs', ...
            'Third argument should be a chebfun2v.');
    end
    
end

if ( nargout > 0 )
    varargout = {h};
end

end

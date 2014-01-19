function varargout = quiver3(Z,F,varargin)
%QUIVER3  3-D quiver plot of a chebfun2v at data mapped by a chebfun2.
%
% QUIVER3(Z,F) plots velocity vectors at the equally spaced surface points
% specified by the chebfun2 Z. We use Z to map a uniform grid. F should be
% a chebfun2v.
%
% QUIVER3(X,Y,Z,F) plots velocity vectors at (x,y,z), where X, Y, Z are
% chebfun2 objects which we use to to map a uniform grid. F should be a
% chebfun2v.
%
% Alternative syntax for this command is:
% QUIVER3(X,Y,Z,[f;g;h]) or QUIVER3(X,Y,Z,f,g,h), where f, g, and h are
% chebfun2 objects.
%
% QUIVER(...,'numpts',N) plots arrows on a N x N uniform grid.
%
%
% This command is a wrapper to CHEBFUN2V/QUIVER3, and is required because
% chebfun2 methods take priority over chebfun2v methods.
%
% See also CHEBFUN2V/QUIVER3.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

numpts = 20; 

if ( ~isa(Z,'chebfun2') )
    error('CHEBFUN2:QUIVER3:INPUT','First argument to this command should be chebfun2.');
end

if ( isempty(varargin) )
    varargin = {};
end

% Number of points to plot
j = 1; argin = {};
while ( ~isempty(varargin) )
    if ( strcmpi(varargin{1}, 'numpts') )
        numpts = varargin{2};
        varargin(1:2) = [];
    else
        argin{j} = varargin{1};
        varargin(1) = [];
        j = j+1;
    end
end
varargin = argin; 



if ( isa(F,'chebfun2v') )              % quiver(Z,F,...)
    zz = new_data_locations(Z, numpts);
    h = quiver3(zz, F, varargin{:});
elseif ( nargin > 3 )
    if ( isa(varargin{2}, 'chebfun2v') ) % quiver(X,Y,Z,F,...)
        xx = new_data_locations(Z, numpts);
        yy = new_data_locations(F, numpts);
        zz = new_data_locations(varargin{1}, numpts);
        h = quiver3( xx, yy, zz, varargin{2}, varargin{3:end} );
    elseif ( isa(Z,'chebfun2') && isa(F,'chebfun2') &&...
            isa(varargin{1},'chebfun2') &&...
            isa(varargin{2},'chebfun2') &&...
            ( nargin ==4 || ~isa(varargin{3},'chebfun2')) )
        FF = vertcat( vertcat(F, varargin{1}), varargin{2}) ;
        h = quiver3( Z, FF );  % call quiver3(Z,F,...)
    elseif ( ( nargin > 5 ) && isa(Z,'chebfun2') && isa(F,'chebfun2') &&...
            isa(varargin{1},'chebfun2') &&...
            isa(varargin{2},'chebfun2') &&...
            isa(varargin{3},'chebfun2') &&...
            isa(varargin{4},'chebfun2') &&...
            ( nargin == 6 || ~isa(varargin{5},'chebfun2')) )
        FF = vertcat( vertcat(varargin{2}, varargin{3}), varargin{4} );
        h = quiver3(Z, F, varargin{1}, FF, varargin{5:end} );  % call quiver3(X,Y,Z,F,...)
    else
        error('CHEBFUN2:QUIVER:INPUTS','Unrecognised input arguments.');
    end
else
    error('CHEBFUN2:QUIVER:INPUTS','Unrecognised input arguments.');
end

if ( nargout > 0 )
    varargout = {h};
end

end


function newloc = new_data_locations(f1, numpts)
% Generate new arrow location if first two inputs are chebfun2 objects.

dom = f1.domain;

% mesh 'em up for the quiver arrows.
x = linspace(dom(1), dom(2), numpts);
y = linspace(dom(3), dom(4), numpts);

[xx, yy] = meshgrid(x, y);
newloc = feval(f1, xx, yy);      % use chebfun2 to generate data locations.

end
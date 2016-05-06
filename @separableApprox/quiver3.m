function varargout = quiver3( Z, F, varargin )
%QUIVER3  3-D quiver plot of a CHEBFUN2V at data mapped by a SEPARABLEAPPROX.
%
% QUIVER3(Z, F) plots velocity vectors at the equally spaced surface points
% specified by the SEPARABLEAPPROX Z. We use Z to map a uniform grid. F should be
% a CHEBFUN2V.
%
% QUIVER3(X, Y, Z, F) plots velocity vectors at (x,y,z), where X, Y, Z are
% SEPARABLEAPPROX objects which we use to to map a uniform grid. F should be a
% CHEBFUN2V.
%
% Alternative syntax for this command is:
% QUIVER3(X,Y,Z,[f;g;h]) or QUIVER3(X,Y,Z,f,g,h), where f, g, and h are
% SEPARABLEAPPROX objects.
%
% QUIVER(...,'numpts',N) plots arrows on a N x N uniform grid.
%
% This command is a wrapper to CHEBFUN2V/QUIVER3, and is required because
% SEPARABLEAPPROX methods take priority over CHEBFUN2V methods.
%
% See also CHEBFUN2V/QUIVER3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

numpts = 20; 

if ( ~isa(Z, 'separableApprox') )
    error('CHEBFUN:SEPARABLEAPPROX:quiver3:inputs1', ...
        'First argument to this command should be SEPARABLEAPPROX.');
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



if ( isa(F, 'chebfun2v') )              % quiver(Z,F,...)
    zz = new_data_locations(Z, numpts);
    h = quiver3(zz, F, varargin{:});
elseif ( nargin > 3 )
    if ( isa(varargin{2}, 'chebfun2v') ) % quiver(X,Y,Z,F,...)
        xx = new_data_locations(Z, numpts);
        yy = new_data_locations(F, numpts);
        zz = new_data_locations(varargin{1}, numpts);
        h = quiver3( xx, yy, zz, varargin{2}, varargin{3:end} );
    elseif ( isa(Z, 'separableApprox') && isa(F, 'separableApprox') &&...
            isa(varargin{1}, 'separableApprox') &&...
            isa(varargin{2}, 'separableApprox') &&...
            ( nargin ==4 || ~isa(varargin{3}, 'separableApprox')) )
        FF = vertcat( vertcat(F, varargin{1}), varargin{2}) ;
        h = quiver3( Z, FF );     % call quiver3(Z,F,...)
    elseif ( ( nargin > 5 ) && isa(Z,'separableApprox') && isa(F,'separableApprox') &&...
            isa(varargin{1}, 'separableApprox') &&...
            isa(varargin{2}, 'separableApprox') &&...
            isa(varargin{3}, 'separableApprox') &&...
            isa(varargin{4}, 'separableApprox') &&...
            ( nargin == 6 || ~isa(varargin{5}, 'separableApprox')) )
        FF = vertcat( vertcat(varargin{2}, varargin{3}), varargin{4} );
        h = quiver3(Z, F, varargin{1}, FF, varargin{5:end} );  % call quiver3(X,Y,Z,F,...)
    else
        error('CHEBFUN:SEPARABLEAPPROX:quiver3:inputs2', ...
            'Unrecognised input arguments.');
    end
else
    error('CHEBFUN:SEPARABLEAPPROX:quiver3:inputs3', 'Unrecognised input arguments.');
end

if ( nargout > 0 )
    varargout = { h };
end

end


function newloc = new_data_locations(f1, numpts)
% Generate new arrow location if first two inputs are SEPARABLEAPPROX objects.

dom = f1.domain;

% mesh 'em up for the quiver arrows.
x = linspace(dom(1), dom(2), numpts);
y = linspace(dom(3), dom(4), numpts);

[xx, yy] = meshgrid(x, y);
newloc = feval(f1, xx, yy);      % use SEPARABLEAPPROX to generate data locations.

end

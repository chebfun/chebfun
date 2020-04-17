function varargout = contour3( f, varargin )
%CONTOUR3   3-D contour plot of a DISKFUN.
%   CONTOUR3(F) is a contour plot of F treating the values of F as heights
%   above the disk. A contour plot shows the level curves of F for some
%   values V. The values V are chosen automatically.
%
%   CONTOUR3(F, N) draws N contour lines, overriding the automatic number.
%   The values V are still chosen automatically.
%   
%   CONTOUR3(F, V) draws LENGTH(V) contour lines at the values specified in
%   the vector V. Use CONTOUR3(F, [V V]) to compute a single contour at the
%   level V.
%
%   CONTOUR3(X, Y, F, ...), CONTOUR3(X, Y, F, N, ...), and
%   CONTOUR3(X, Y, F, V, ...) use matrices X and Y to specify the plotting
%   grid.
%
%   [C, H] = CONTOUR3(...) returns contour matrix C as described in
%   CONTOURC and a handle H to a contour object. This handle can be used as
%   input to CLABEL.
%
%   CONTOUR3(F, 'NUMPTS', N) plots the contour lines on an N by N uniform
%   grid. If NUMPTS is not given then we plot on a 200 by 200 grid.
%
%   CONTOUR3(F, 'PIVOTS', STR) plots the contour lines with the pivot
%   locations used during construction.
%
% See also CONTOUR, CONTOURF.

% Copyright 2020 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty( f ) )  % Empty check.
    contour3( [] );
    return
end

% Minimum number of plotting points.
minplotnum = 200;
doPivotPlot = 0; 

% Extract from the inputs the user defined options.
j = 1; 
argin = {};
while ( ~isempty( varargin ) )
    if ( strcmpi( varargin{1}, 'numpts') ) % If given numpts then use them.
        minplotnum = varargin{2};
        varargin(1:2) = [];
    elseif ( strcmpi(varargin{1}, 'pivots') ) % Should we plot the pivots?
        doPivotPlot = 1;
        if ( length( varargin ) < 2 ) 
            error('CHEBFUN:DISKFUN:contour3:pivotStyle', ...
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
numargs = numel(argin);

% Did the user want a plot of the pivot locations?
if ( doPivotPlot )
    if ( ( ~isempty(argin) ) && ( length(argin{1}) < 5 ) )
        % Column, row, pivot plot
        holdState = ishold;
        plot( f, argin{:} ), hold on
        argin(1) = [];
        contour3( f, argin{:} )
        if ( ~holdState )
            hold off
        end
        return
    end
end

if ( isa(f, 'double') )
    % CONTOUR3(xx, yy, F,...)

    if ( (numargs >= 2) && isa(argin{1}, 'double') && isa(argin{2}, 'diskfun') )
        % Extract inputs:
        xx = f; 
        yy = argin{1}; 
        f = argin{2};
        argin(1:2) = [];
        % Evaluate diskfun: 
        vals = feval(f, xx, yy, 'polar');

    else
        error('CHEBFUN:DISKFUN:contour3:badInputs', ...
            'Unrecognised input arguments.');
    end

elseif ( isa(f, 'diskfun') )

    dom = f.domain;
    if ( (numargs >= 2) && isa(argin{1}, 'diskfun') && isa(argin{2}, 'diskfun') )
        % CONTOUR3 plot on a surface.

        % Extract inputs:
        xx = f; 
        yy = argin{1}; 
        f = argin{2};
        argin(1:2) = [];

        % Check CONTOUR3 objects are on the same domain.
        if ( ~domainCheck(xx, yy) )
            error('CHEBFUN:DISKFUN:contour3:domains', ...
                'Domains of DISKFUN objects are not consistent.');
        end

        % Evaluate f on equally spaced grid:
        x = linspace( dom(1), dom(2), minplotnum );
        y = linspace( dom(3), dom(4), minplotnum );
        [mxx, myy] = meshgrid(x, y);
        xx = feval(xx, mxx, myy, 'polar'); 
        yy = feval(yy, mxx, myy, 'polar');
        [xx,yy] = cart2pol(xx,yy);
        vals = feval(f, xx, yy, 'polar');

    elseif ( (numargs == 2) || ((numargs > 2) && ~isa(argin{1}, 'diskfun')) ) 
        % CONTOUR3(xx, yy, f)

        % Evaluate f at equally spaced points.
        x = linspace( dom(1), dom(2), minplotnum );
        y = linspace( dom(3), dom(4), minplotnum );
        [xx, yy] = meshgrid(x, y);
        vals = feval(f, xx, yy, 'polar');

    elseif ( (numargs == 0) || ((numargs > 0) && isa(argin{1}, 'double')) )
        % CONTOUR3(f) 

        % Evaluate at equally spaced grid: 
        x = linspace( dom(1), dom(2), minplotnum );
        y = linspace( dom(3), dom(4), minplotnum );
        [xx, yy] = meshgrid(x, y);
        vals = feval(f, xx, yy, 'polar');

    else
        error('CHEBFUN:DISKFUN:contour3:inputs1', ...
            'Unrecognised input arguments.');
    end

else

    error('CHEBFUN:DISKFUN:contour3:inputs2', ...
        'Unrecognised input arguments.');

end

[X, Y] = pol2cart(xx, yy);
% CONTOUR3 plot:
[c, h] = contour3( X, Y, vals, argin{:} );
axis square

% Return plot handle if appropriate.
if ( nargout == 1 )
    varargout = {h};
elseif ( nargout == 2 )
    varargout = {c, h};
end

end

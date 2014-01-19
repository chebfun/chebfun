function varargout = contourf( f, varargin )
%CONTOURF Filled contour plot of a chebfun2.
%
% CONTOURF(...) is the same as CONTOUR(...) except that the areas
% between contours are filled with colors according to the Z-value
% for each level.  Contour regions with data values at or above a
% given level are filled with the color that maps to the interval.
%
% NaN's in the Z-data leave white holes with black borders in the
% contour plot.
%
% When you use the CONTOURF(Z, V) syntax to specify a vector of contour
% levels (V must increase monotonically), contour regions with
% Z-values less than V(1) are not filled (are rendered in white).
% To fill such regions with a color, make V(1) less than or equal to
% the minimum Z-data value.
%
% CONTOURF(F,'NUMPTS',N) computes the contour lines on a N by N grid. If N
%     is larger than 200 then the contour lines are drawn with more detail.
%
% [C, H] = CONTOURF(...) also returns a handle H to a CONTOURGROUP object.
%
% See also CONTOUR.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty(f) )  % empty check.
    contourf([]);
    return
end

% Min plotting number of points:
minplotnum = 200;
dom = f.domain;
doPivotPlot = 0; 

% Number of points to plot
j = 1; argin = {}; pivots = 0;
while ( ~isempty(varargin) )
    if ( strcmpi(varargin{1},'numpts') ) % If given numpts then use them.
        minplotnum = varargin{2};
        varargin(1:2) = [];
    elseif ( strcmpi(varargin{1},'pivots') ) % If given numpts then use them.
        doPivotPlot = 1;
        argin{j} = varargin{2};
        varargin(1:2) = [];
    else
        argin{j} = varargin{1};
        varargin(1) = [];
        j = j+1;
    end
end

if ( doPivotPlot )    % Do pivot plot. 
    if ( ( ~isempty(argin) ) && ( length(argin{1}) < 5 ) )
        % Column, row, pivot plot
        plot( f, argin{:} ), hold on
        argin(1) = [];
        contourf( f, argin{:} ), hold off
        return
    end
end


if ( isa(f,'double') )   
    % contourf(xx,yy,F,...)
    if ( ( nargin >= 3 ) && ( isa(varargin{1},'double') ) && ...
                                                 isa(varargin{2},'chebfun2') )
        % Extract inputs:
        xx = f; 
        yy = varargin{1}; 
        f = varargin{2};
        % Evaluate chebfun2: 
        vals = feval(f, xx, yy);
        % contourf plot:
        [c, h] = contourf( xx, yy, vals, varargin{3:end} );
    else
        error('CHEBFUN2:contourf:INPUTS','Unrecognised input arguments.');
    end
    
elseif ( isa(f,'chebfun2') )     
    % contourf( f )
    if ( ( nargin == 3 ) || ( ( nargin > 3 ) && ( ~isa(argin{1},'chebfun2') ) ) )
        % Evaluate f at equally spaced points.
        x = linspace( dom(1), dom(2), minplotnum );
        y = linspace( dom(3), dom(4), minplotnum );
        [xx, yy] = meshgrid(x, y);
        vals = feval( f, xx, yy );
        % contourf plot:
        [c, h] = contourf( xx, yy, vals, argin{:} );
    elseif ( ( nargin >= 3 ) && ( isa(argin{1},'chebfun2') ) && ...  % contourf plot on surface.
                                                  ( isa(argin{2},'chebfun2') ) )
        % Extract inputs:
        xx = f; 
        yy = argin{1}; 
        f = argin{2};
        % check chebfun2 objects are on the same domain.
        dom = xx.domain; domcheck = yy.domain; domf = f.domain;
        if any(dom - domcheck) || any(dom-domf)
            error('CHEBFUN2:contourf:DOMAIN','Domains of chebfun2 objects are not consistent.');
        end
        % Evaluate f on equally spaced grid:
        x = linspace( dom(1), dom(2), minplotnum );
        y = linspace( dom(3), dom(4), minplotnum );
        [mxx, myy] = meshgrid(x, y);
        xx = feval(xx, mxx, myy); 
        yy = feval(yy, mxx, myy);
        vals = feval(f, mxx, myy);
        % contourf plot:
        [c,h] = contourf( xx, yy, vals, argin{3:end} );
    elseif ( ( nargin == 1) || ( ( nargin > 1 ) && ( isa(argin{1},'double') ) ) )
        % Evaluate at equally spaced grid: 
        x = linspace( dom(1), dom(2), minplotnum );
        y = linspace( dom(3), dom(4), minplotnum );
        [xx, yy] = meshgrid(x, y);
        vals = feval(f, xx, yy );
        % contourf plot:
        [c, h] = contourf( xx, yy, vals, argin{:} );
    else
        error('CHEBFUN2:contourf:INPUTS','Unrecognised input arguments.');
    end
else
    error('CHEBFUN2:contourf:INPUTS','Unrecognised input arguments.');
end

% return plot handle if appropriate.
if ( nargout > 0 )
    varargout = {h};
end
if ( nargout > 1 )
    varargout = {c,h};
end
end
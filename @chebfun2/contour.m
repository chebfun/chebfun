function varargout = contour( f, varargin )
%CONTOUR  contour plot of a chebfun2.
%
% CONTOUR(F) is a contour plot of F treating the values of F as heights
% above a plane. A contour plot are the level curves of F for some values
% V. The values V are chosen automatically.
%
% CONTOUR(F,N) draw N contour lines, overriding the automatic number. The
% values V are still chosen automatically.
%
% CONTOUR(F, V) draw LENGTH(V) contour lines at the values specified in
% vector V.  Use contour(F, [v, v]) to compute a single contour at the
% level v.
%
% CONTOUR(X,Y,F,...), CONTOUR(X,Y,F,N,...), and CONTOUR(X,Y,F,V,...) where
% X and Y are matrices or chebfun2 objects that are used to specify the
% plotting grid.
%
% [C, H] = contour(...) returns contour matrix C as described in
% CONTOURC and a handle H to a contourgroup object.  This handle can
% be used as input to CLABEL.
%
% CONTOUR(F,'NUMPTS',N) plots the contour lines on a N by N uniform grid.
% If NUMPTS is not given then we plot at a N by N grid, where N is decided
% by a field in the chebfun2pref object.
%
% CONTOUR(F,'PIVOTS',STR) plots the contour lines with the pivot locations used
% during constructor.
%
% See also CONTOURF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty(f) )  % empty check.
    contour([]);
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
        contour( f, argin{:} ), hold off
        return
    end
end


if ( isa(f,'double') )   
    % contour(xx,yy,F,...)
    if ( ( nargin >= 3 ) && ( isa(varargin{1},'double') ) && ...
                                                 isa(varargin{2},'chebfun2') )
        % Extract inputs:
        xx = f; 
        yy = varargin{1}; 
        f = varargin{2};
        % Evaluate chebfun2: 
        vals = feval(f, xx, yy);
        % contour plot:
        [c, h] = contour( xx, yy, vals, varargin{3:end} );
    else
        error('CHEBFUN2:CONTOUR:INPUTS','Unrecognised input arguments.');
    end
    
elseif ( isa(f,'chebfun2') )     
    % contour( f )
    if ( ( nargin == 3 ) || ( ( nargin > 3 ) && ( ~isa(argin{1},'chebfun2') ) ) )
        % Evaluate f at equally spaced points.
        x = linspace( dom(1), dom(2), minplotnum );
        y = linspace( dom(3), dom(4), minplotnum );
        [xx, yy] = meshgrid(x, y);
        vals = feval( f, xx, yy );
        % contour plot:
        [c, h] = contour( xx, yy, vals, argin{:} );
    elseif ( ( nargin >= 3 ) && ( isa(argin{1},'chebfun2') ) && ...  % contour plot on surface.
                                                  ( isa(argin{2},'chebfun2') ) )
        % Extract inputs:
        xx = f; 
        yy = argin{1}; 
        f = argin{2};
        % check chebfun2 objects are on the same domain.
        dom = xx.domain; domcheck = yy.domain; domf = f.domain;
        if any(dom - domcheck) || any(dom-domf)
            error('CHEBFUN2:CONTOUR:DOMAIN','Domains of chebfun2 objects are not consistent.');
        end
        % Evaluate f on equally spaced grid:
        x = linspace( dom(1), dom(2), minplotnum );
        y = linspace( dom(3), dom(4), minplotnum );
        [mxx, myy] = meshgrid(x, y);
        xx = feval(xx, mxx, myy); 
        yy = feval(yy, mxx, myy);
        vals = feval(f, mxx, myy);
        % contour plot:
        [c,h] = contour( xx, yy, vals, argin{3:end} );
    elseif ( ( nargin == 1) || ( ( nargin > 1 ) && ( isa(argin{1},'double') ) ) )
        % Evaluate at equally spaced grid: 
        x = linspace( dom(1), dom(2), minplotnum );
        y = linspace( dom(3), dom(4), minplotnum );
        [xx, yy] = meshgrid(x, y);
        vals = feval(f, xx, yy );
        % contour plot:
        [c, h] = contour( xx, yy, vals, argin{:} );
    else
        error('CHEBFUN2:CONTOUR:INPUTS','Unrecognised input arguments.');
    end
else
    error('CHEBFUN2:CONTOUR:INPUTS','Unrecognised input arguments.');
end

% return plot handle if appropriate.
if ( nargout > 0 )
    varargout = {h};
end
if ( nargout > 1 )
    varargout = {c,h};
end
end
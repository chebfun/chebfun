function varargout = contourf( f, varargin )
%CONTOURF   Filled contour plot of a SEPARABLEAPPROX.
%   CONTOURF(...) is the same as CONTOUR(...) except that the areas between
%   contours are filled with colors according to the Z-value for each level.
%   Contour regions with data values at or above a given level are filled with
%   the color that maps to the interval.
%
%   NaN's in the Z-data leave white holes with black borders in the contour
%   plot.
%
%   When you use the CONTOURF(Z, V) syntax to specify a vector of contourf
%   levels (V must increase monotonically), contourf regions with Z-values less
%   than V(1) are not filled (are rendered in white). To fill such regions with
%   a color, make V(1) less than or equal to the minimum Z-data value.
%
%   CONTOURF(F, 'NUMPTS', N) computes the contourf lines on a N by N grid. If N
%   is larger than 200 then the contourf lines are drawn with more detail.
%
%   [C, H] = CONTOURF(...) also returns a handle H to a CONTOURGROUP object.
%
% See also CONTOUR.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Refactor to reduce duplication.

if ( isempty( f ) )  % Empty check.
    contourf( [] );
    return
end

% Minimum number of plotting points:
minplotnum = 200;
doPivotPlot = 0; 

% Extract from the inputs the user defined options: 
j = 1; 
argin = {};
while ( ~isempty( varargin ) )
    if ( strcmpi( varargin{1}, 'numpts') ) % If given numpts then use them.
        minplotnum = varargin{2};
        varargin(1:2) = [];
    elseif ( strcmpi(varargin{1}, 'pivots') ) % If given numpts then use them.
        doPivotPlot = 1;
        if ( length( varargin ) < 2 ) 
            error('CHEBFUN:SEPARABLEAPPROX:contourf:pivotStyle', ...
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
    if ( ( ~isempty(argin) ) && ( length(argin{1}) < 5 ) )
        % Column, row, pivot plot
        plot( f, argin{:} ), hold on
        argin(1) = [];
        contourf( f, argin{:} ), hold off
        return
    end
end

if ( isa(f, 'double') )                
    % CONTOUR(xx, yy, F,...)
    
    if ( (nargin >= 3) && isa(argin{1}, 'double') && isa(argin{2}, 'separableApprox') )
        % Extract inputs:
        xx = f; 
        yy = argin{1}; 
        f = argin{2};
        argin(1:2) = [];
        % Evaluate separableApprox: 
        vals = feval(f, xx, yy);
        
    else
        error('CHEBFUN:SEPARABLEAPPROX:contourf:inputs1', ...
            'Unrecognised input arguments.');
    end
    
elseif ( isa(f, 'separableApprox') ) 
    
    dom = f.domain;
    if ( (nargin == 3) || (nargin > 3) && ~isa(argin{1},'separableApprox') ) 
        % CONTOUR(xx, yy, f)
        
        % Evaluate f at equally spaced points.
        x = linspace( dom(1), dom(2), minplotnum );
        y = linspace( dom(3), dom(4), minplotnum );
        [xx, yy] = meshgrid(x, y);
        vals = feval( f, xx, yy );

    elseif ( (nargin >= 3) && isa(argin{1},'separableApprox') && isa(argin{2},'separableApprox') )
        % CONTOUR plot on a surface.
        
        % Extract inputs:
        xx = f; 
        yy = argin{1}; 
        f = argin{2};
        argin(1:2) = [];
        
        % Check CONTOUR objects are on the same domain.
        if ( ~domainCheck(xx, yy) || ~domainCheck(yy, f) )
            error('CHEBFUN:SEPARABLEAPPROX:contourf:domainMismatch', ...
                'Domains of SEPARABLEAPPROX objects are not consistent.');
        end
        
        % Evaluate f on equally spaced grid:
        x = linspace( dom(1), dom(2), minplotnum );
        y = linspace( dom(3), dom(4), minplotnum );
        [mxx, myy] = meshgrid(x, y);
        xx = feval(xx, mxx, myy); 
        yy = feval(yy, mxx, myy);
        vals = feval(f, mxx, myy);

    elseif ( ( nargin == 1) || ( ( nargin > 1 ) && ( isa(argin{1},'double') ) ) )    
        % CONTOUR(f) 
        
        % Evaluate at equally spaced grid: 
        x = linspace( dom(1), dom(2), minplotnum );
        y = linspace( dom(3), dom(4), minplotnum );
        [xx, yy] = meshgrid(x, y);
        vals = feval(f, xx, yy );
        
    else
        error('CHEBFUN:SEPARABLEAPPROX:contourf:inputs2', ...
            'Unrecognised input arguments.');
    end
    
else
    
    error('CHEBFUN:SEPARABLEAPPROX:contourf:inputs3', ...
        'Unrecognised input arguments.');
    
end

% Contour plot:
[c, h] = contourf( xx, yy, vals, argin{:} );

% return plot handle if appropriate.
if ( nargout == 1 )
    varargout = {h};
elseif ( nargout == 2 )
    varargout = {c, h};
end

end

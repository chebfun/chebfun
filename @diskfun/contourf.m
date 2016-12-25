function varargout = contourf( f, varargin )
%CONTOURF   Filled contour plot of a DISKFUN.
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

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


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
            error('CHEBFUN:DISKFUN:contourf:pivotStyle', ...
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
    error('CHEBFUN:DISKFUN:contourf:pivotstyle', ...
            'Pivots cannot be plotted with ''contourf''. Use ''contour'' instead.');
end

if ( isa(f, 'double') )                
    % CONTOURF(xx, yy, F,...)
    
    if ( (nargin >= 3) && isa(argin{1}, 'double') && isa(argin{2}, 'diskfun') )
        % Extract inputs:
        xx = f; 
        yy = argin{1}; 
        f = argin{2};
        argin(1:2) = [];
        % Evaluate separableApprox: 
        vals = feval(f, xx, yy, 'polar');
        
    else
        error('CHEBFUN:DISKFUN:contourf:badInputs', ...
            'Unrecognised input arguments.');
    end
    
elseif ( isa(f, 'diskfun') ) 
    
    dom = f.domain;
    if ( (nargin >= 3) && isa(argin{1},'diskfun') && isa(argin{2},'diskfun') )
        % CONTOURF plot on a surface.
        
        % Extract inputs:
        xx = f; 
        yy = argin{1}; 
        f = argin{2};
        argin(1:2) = [];
        
        % Check CONTOUR objects are on the same domain.
        if ( ~domainCheck(xx, yy) )
            error('CHEBFUN:DISKFUN:contourf:domains', ...
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
        
    elseif ( (nargin == 3) || (nargin > 3) && ~isa(argin{1},'diskfun') ) 
        % CONTOURF(xx, yy, f)
        
        % Evaluate f at equally spaced points.
        x = linspace( dom(1), dom(2), minplotnum );
        y = linspace( dom(3), dom(4), minplotnum );
        [xx, yy] = meshgrid(x, y);
        vals = feval( f, xx, yy, 'polar' );
        
    elseif ( ( nargin == 1) || ( ( nargin > 1 ) && ( isa(argin{1},'double') ) ) )    
        % CONTOURF(f) 
        
        % Evaluate at equally spaced grid: 
        x = linspace( dom(1), dom(2), minplotnum );
        y = linspace( dom(3), dom(4), minplotnum );
        [xx, yy] = meshgrid(x, y);
        vals = feval(f, xx, yy, 'polar' );
        % This code would be faster:
%         % Evaluate at Fourier-Chebyshev grid: 
%         x = [trigpts( minplotnum-1, [dom(1) dom(2)] ); dom(1)];
%         y = chebpts(2*minplotnum-1,[-1 1]); y = y(minplotnum:end);
%         [xx, yy] = meshgrid(x, y);
%         vals = sample(f, minplotnum-1, minplotnum);
%         vals = [vals vals(:,1)];
    else
        error('CHEBFUN:DISKFUN:contourf:inputs1', ...
            'Unrecognised input arguments.');
    end
else
    error('CHEBFUN:DISKFUN:contourf:inputs2', ...
        'Unrecognised input arguments.');
end

[X, Y] = pol2cart(xx, yy); 
% Contour plot:
[c, h] = contourf( X, Y, vals, argin{:} );
axis square

% Add unit circle
holdState = ishold; 
circ = exp(1i*pi*linspace(-1,1,101));
hold on 
plot(real(circ), imag(circ), 'k-', 'Linewidth', .3)
hold off

if holdState
    hold on;
else
    hold off;
end

% Return plot handle if appropriate.
if ( nargout == 1 )
    varargout = {h};
elseif ( nargout == 2 )
    varargout = {c, h};
end

end

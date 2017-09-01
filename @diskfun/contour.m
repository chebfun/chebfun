function varargout = contour( f, varargin )
%CONTOUR   contour plot of a DISKFUN.
%   CONTOUR(F) is a contour plot of F treating the values of F as heights
%   above the disk. A contour plot shows the level curves of F for some
%   values V. The values V are chosen automatically.
%
%   CONTOUR(F, N) draw N contour lines, overriding the automatic number.
%   The values V are still chosen automatically.
%   
%   CONTOUR(F, V) draws LENGTH(V) contour lines at the values specified in
%   vector V. Use contour(F, [v, v]) to compute a single contour at the
%   level v.
%  
%   [C, H] = contour(...) returns contour matrix C as described in CONTOURC
%   and a handle H to a contourgroup object.  This handle can be used as
%   input to CLABEL.
%
%   CONTOUR(F, 'NUMPTS', N) plots the contour lines on a N by N uniform
%   grid. If NUMPTS is not given then we plot on an 200 by 200 grid.
%
%   CONTOUR(F, 'PIVOTS', STR) plots the contour lines with the pivot locations
%   used during constructor.
%
% See also CONTOURF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty( f ) )  % Empty check.
    contour( [] );
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
    elseif ( strcmpi(varargin{1}, 'pivots') ) % If the pivots are to be plotted
        doPivotPlot = 1;
        if ( length( varargin ) < 2 ) 
            error('CHEBFUN:DISKFUN:contour:pivotStyle', ...
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
if ( doPivotPlot )    % Do pivot plot. 
    if ( ( ~isempty(argin) ) && ( length(argin{1}) < 5 ) )
        % Column, row, pivot plot
        holdState = ishold;
        plot( f, argin{:} ), hold on
        argin(1) = [];
        contour( f, argin{:} )
        if ( ~holdState )
            hold off;
        end
        return
    end
end

if ( isa(f, 'double') )                
    % CONTOUR(xx, yy, F,...)
    
    if ( (numargs >= 2) && isa(argin{1}, 'double') && isa(argin{2}, 'diskfun') )
        % Extract inputs:
        xx = f; 
        yy = argin{1}; 
        f = argin{2};
        argin(1:2) = [];
        % Evaluate separableApprox: 
        vals = feval(f, xx, yy, 'polar');
        
    else
        error('CHEBFUN:DISKFUN:contour:badInputs', ...
            'Unrecognised input arguments.');
    end
    
elseif ( isa(f, 'diskfun') ) 
    
    dom = f.domain;
    if ( (numargs >= 2) && isa(argin{1},'diskfun') && isa(argin{2},'diskfun') )
        % CONTOUR plot on a surface.
        
        % Extract inputs:
        xx = f; 
        yy = argin{1}; 
        f = argin{2};
        argin(1:2) = [];
        
        % Check CONTOUR objects are on the same domain.
        if ( ~domainCheck(xx, yy) )
            error('CHEBFUN:DISKFUN:contour:domains', ...
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
        
    elseif ( (numargs == 2) || (numargs > 2) && ~isa(argin{1},'diskfun') ) 
        % CONTOUR(xx, yy, f)
        
        % Evaluate f at equally spaced points.
        x = linspace( dom(1), dom(2), minplotnum );
        y = linspace( dom(3), dom(4), minplotnum );
        [xx, yy] = meshgrid(x, y);
        vals = feval( f, xx, yy, 'polar' );
        
    elseif ( ( numargs == 0) || ( ( numargs > 0 ) && ( isa(argin{1},'double') ) ) )    
        % CONTOUR(f) 
        
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
        error('CHEBFUN:DISKFUN:contour:inputs1', ...
            'Unrecognised input arguments.');
    end
else
    error('CHEBFUN:DISKFUN:contour:inputs2', ...
        'Unrecognised input arguments.');
end

[X, Y] = pol2cart(xx, yy); 
% Contour plot:
[c, h] = contour( X, Y, vals, argin{:} );
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

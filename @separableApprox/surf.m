function varargout = surf( f, varargin )
%SURF  Surface plot of a SEPARABLEAPPROX.
%   SURF(F, C) plots the colored parametric surface defined by F and the matrix
%   C. The matrix C, defines the colouring of the surface.
%
%   SURF(F) uses colors proportional to surface height.
%
%   SURF(X, Y, F, ...) is the same as SURF(F, ...) when X and Y are SEPARABLEAPPROX
%   objects except X and Y supplies the plotting locations are  mapped by X and
%   Y.
%
%   SURF(..., 'PropertyName', PropertyValue,...) sets the value of the specified
%   surface property. Multiple property values can be set with a single
%   statement.
%
%   H = SURF(...) returns a handle to a surface plot object.
%
% See also PLOT, SURFC.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(f) )
    h = surf([]);
    if ( nargout == 1 )
        varargout = {h};
    end
    return
end

% How dense to make the samples.
minPlotNum = 200;
defaultOpts = {'facecolor', 'interp', 'edgealpha', .5, 'edgecolor', 'none'};

% Number of points to plot
j = 1; argin = {};
while ( ~isempty(varargin) )
    if strcmpi(varargin{1}, 'numpts')
        minPlotNum = varargin{2};
        varargin(1:2) = [];
    else
        argin{j} = varargin{1};
        varargin(1) = [];
        j = j+1;
    end
end

if ( isempty(argin) )
    argin = {};
end

if ( isa(f,'separableApprox') )
    
    if ( (nargin == 1) || ...
            ( (nargin > 1) && ~isempty(argin) && ~isa(argin{1}, 'separableApprox') ) || ...
            ( (nargin == 3) && isempty(argin)) )
        % surf(f,...)
        
        dom = f.domain;
        xdata = linspace(dom(1), dom(2), minPlotNum);
        ydata = linspace(dom(3), dom(4), minPlotNum);
        [xx, yy] = meshgrid(xdata, ydata);
        C = feval(f, xx, yy); 
        % Make some corrections to C for prettier plotting.
        if ( norm(C - C(1,1),inf) < 1e-10 )
            % If vals are very close up to round off then the color scale is
            % hugely distorted. This fixes that.
            [n, m] = size(C);
            C = C(1,1)*ones(n, m);
        end
        h = surf(xx, yy, C, defaultOpts{:}, argin{:});
        xlim(dom(1:2)), ylim(dom(3:4))
        
    elseif ( nargin > 2)
        % surf(x,y,f,...), with x, y, f SEPARABLEAPPROX object
        
        x = f;
        y = argin{1};
        if ( isa(y, 'separableApprox') )
            % Check domains of x and y are the same.
            dom = x.domain;
            rectcheck = y.domain;
            if ( any(dom - rectcheck) )
                error('CHEBFUN:SEPARABLEAPPROX:surf:domainMismatch', ...
                    'Domains of SEPARABLEAPPROX objects do not match.');
            end
        end
        xdata = linspace(dom(1), dom(2), minPlotNum);
        ydata = linspace(dom(3), dom(4), minPlotNum);
        [xx, yy] = meshgrid(xdata, ydata);
        x = feval(x, xx, yy);
        y = feval(y, xx, yy);
        if ( isa(argin{2}, 'separableApprox') )         % surf(x,y,f,...)
            vals = feval(argin{2}, xx, yy);
            if ( nargin < 4 )                    % surf(x,y,f)
                C = vals;
            elseif ( isa(argin{3}, 'double') )   % surf(x,y,f,C,...)
                C = argin{3};
                argin(3) = [];
            elseif ( isa(argin{3}, 'separableApprox'))  % Colour matrix is a SEPARABLEAPPROX.
                C = feval(argin{3},xx,yy);
                argin(3) = [];
            else
                C = vals;
            end
            
            % Make some corrections to C for prettier plotting.
            if ( norm(C - C(1,1),inf) < 1e-10 )
                % If vals are very close up to round off then the color scale is
                % hugely distorted. This fixes that.
                [n, m] = size(C);
                C = C(1,1)*ones(n, m);
            end
            
            h = surf(x, y, vals, C, defaultOpts{:}, argin{3:end});
            xlabel('x'), ylabel('y')
            
            % There is a bug in matlab surf plot when vals are very nearly a
            % constant. Fix this manually by resetting axis scaling.
            if ( norm(vals - vals(1,1),inf) < 1e-10*norm(vals,inf) && ...
                    norm(vals - vals(1,1),inf) > 0 )
                v = vals(1,1);
                absv = abs(v);
                zlim([v-.5*absv v+.5*absv])
            end
        else
            error('CHEBFUN:SEPARABLEAPPROX:surf:inputs1', ...
                ['The third argument should be a SEPARABLEAPPROX', ...
                'if you want to supply SEPARABLEAPPROX data.'])
        end
        
    else  %surf(f,C)
        dom = f.domain;
        x = chebfun2(@(x,y) x,dom);
        y = chebfun2(@(x,y) y,dom);
        h = surf(x, y, f, argin{1}, defaultOpts{:}, argin{2:end});
        xlim(dom(1:2)), ylim(dom(3:4))
        
    end
    
else     % surf(X,Y,f,...)
    error('CHEBFUN:SEPARABLEAPPROX:surf:inputs2', ...
        ['Data should be given as SEPARABLEAPPROX objects \n ', ...
        'For example: \n x = chebfun2(@(x,y)x); ', ...
        'y = chebfun2(@(x,y)y);\n surf(x,y,f)']);
    
end

if ( nargout > 0 )
    varargout = {h};
end

end

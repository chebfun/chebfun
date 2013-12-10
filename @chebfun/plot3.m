function varargout = plot3(f, varargin)
%PLOT3   Plot for CHEBFUN objects in 3-D space. 
%   PLOT3() is a three-dimensional analogue of PLOT().
%   
%   PLOT3(X, Y, Z), where X, Y, and Z are three CHEBFUN objects, plots a line in
%   3-space. X, Y, and Z may be array-valued, but must have the same number of
%   columns.
%
%   PLOT3(F), where F is an array-valued CHEBFUN, plots a 3D plot where the kth
%   column of F is plotted along the line y = k.
%   
%   Various line types, plot symbols, and colors may be obtained with PLOT3(X,
%   Y, Z, S) or PLOT3(F), where S is a string of length 1, 2 or 3 containing
%   characters listed under the PLOT command.
%
%   [HLINE, HPOINT, HJUMP] = PLOT(X, Y, Z) returns column vectors of handles to
%   lineseries objects (one for each column in the case of array-valued CHEBFUN
%   objects), corresponding to the line, point, and jump plots, respectively.
%
% Example: A helix:      
%   t = chebfun('t', [0, 10*pi]);
%   plot3(sin(t), cos(t), t);
% Example: Chebyshev polynomials:
%   T = chebpoly(0:5);
%   plot3(T);
%
% See also PLOT, PLOTDATA.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% [TODO]: Implement plotting of delta functions.

% Deal with an empty input:
if ( isempty(f) )
    if ( nargout == 1 )
        varargout{1} = plot3([]);
    end
    return
end

% Store the hold state of the current axis:
holdState = ishold;

%%

if ( nargin == 1 || ~isa(varargin{2}, 'chebfun') )
    % PLOT(F): (Plot an array-valued CHEBFUN in 3D.)
    
    % TODO: This should be moved to WATERFALL().
    
    % Get the data from PLOTDATA():
    data = plotData(f);
    numCols = min(size(f));
    
    % Duplicate the x data for each column:
    data.xLine = repmat(data.xLine, 1, numCols);
    data.xPoints = repmat(data.xPoints, 1, numCols);
    data.xJumps = repmat(data.xJumps, 1, numCols);
    % Copy the y data to z:
    data.zLine = data.yLine;
    data.zPoints = data.yPoints;
    data.zJumps = data.yJumps;
    % Make discrete data in the y direction:
    data.yLine = repmat(1:numCols, size(data.xLine, 1), 1);
    data.yPoints = repmat(1:numCols, size(data.xPoints, 1), 1);
    data.yJumps = repmat(1:numCols, size(data.xJumps, 1), 1);

else
    % PLOT(F, G, H):
    
    % Extract the other CHEBFUN objects from varargin:
    g = varargin{1};
    h = varargin{2};
    varargin(1:2) = [];
    
    % We can only plot real against real:
    if ( ~isreal(f) || ~isreal(g) || ~isreal(h) )
        warning('CHEBFUN:plot:complex', ...
            'Warning: Imaginary parts of complex X, Y, and/or Z arguments ignored.');
    end
    f = real(f);
    g = real(g);
    h = real(h);
    
    % Get the data for plotting from PLOTDATA():
    data = plotData(f, g, h);
end


%%
% Plot the curve
h1 = plot3(data.xLine, data.yLine, data.zLine, varargin{:});
set(h1, 'Marker', 'none')
hold on

%%
% Plot the points:
h2 = plot3(data.xPoints, data.yPoints, data.zPoints, varargin{:});

% Change the style accordingly:
set(h2, 'LineStyle', 'none')

%%
% Plot the jumps:;
h3 = plot3(data.xJumps, data.yJumps, data.zJumps, varargin{:});

% Change the style accordingly:
set(h3,'LineStyle', ':', 'Marker', 'none')

%%
% Return hold state to what it was before:
if ( ~holdState )
    hold off
end

% Give an output if one was requested:
if ( nargout > 0 )
    varargout{1} = h1;
    varargout{2} = h2;
    varargout{3} = h3;
end

end
function varargout = plot3(f, g, h, varargin)
%PLOT3   Plot for CHEBFUN objects in 3-D space.
%   PLOT3 is a three-dimensional analogue of PLOT.
%
%   PLOT3(X, Y, Z), where X, Y, and Z are three CHEBFUN objects, plots a line in
%   3-space. X, Y, and Z may be array-valued, but must have the same number of
%   columns.
%
%   Various line types, plot symbols, and colors may be obtained with PLOT3(X,
%   Y, Z, S) where S is a string of length 1, 2 or 3 containing characters
%   listed under the PLOT command.
%
%   [HLINE, HPOINT, HJUMP] = PLOT(X, Y, Z) returns column vectors of handles to
%   lineseries objects (one for each column in the case of array-valued CHEBFUN
%   objects), corresponding to the line, point, and jump plots, respectively.
%
% Example: A helix:
%   t = chebfun('t', [0, 10*pi]);
%   plot3(sin(t), cos(t), t);
%
% See also PLOT, PLOTDATA.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% Deal with an empty input:
if ( isempty(f) || isempty(g) || isempty(h) )
    if ( nargout == 1 )
        varargout{1} = plot3([]);
    end
    return
end

%% Initialization:

% Suppress inevitable warning for growing these arrays:
%#ok<*AGROW>

% Store the hold state of the current axis:
holdState = ishold;

% We can only plot real against real:
if ( ~isreal(f) || ~isreal(g) || ~isreal(h) )
    warning('CHEBFUN:CHEBFUN:plot:complex', ...
        'Warning: Imaginary parts of complex X, Y, and/or Z arguments ignored.');
end
f = real(f);
g = real(g);
h = real(h);

% Call PLOTDATA():
if ( numel(f) == 1 && numel(g) == 1 && numel(h) == 1)
    % Array-valued CHEBFUN case:
    newData = plotData(f, g, h);
else
    % QUASIMATRIX case:
    numCols = max([numColumns(f), numColumns(g), numColumns(h)]);
    f = expand(num2cell(f), numCols);
    g = expand(num2cell(g), numCols);
    h = expand(num2cell(h), numCols);
    for k = numel(f):-1:1
        newData(k) = plotData(f{k}, g{k}, h{k});
    end
end

% Remove global plotting options from input arguments.
[lineStyle, pointStyle, jumpStyle, deltaStyle, varargin] = ...
    chebfun.parsePlotStyle(varargin{:});

% Style data.
pos = 0; styleData = [];
% Find the location of the next CHEBFUN in the input array:
while ( (pos < length(varargin)) && ~isa(varargin{pos + 1}, 'chebfun') )
    pos = pos + 1;
end
if ( pos > 0 )
    styleData = varargin(1:pos);
    varargin(1:pos) = [];
end

% Initialise storage:
lineData = {};
pointData = {};
jumpData = {};

%% Preparation of the data:

% Loop over the columns:
for k = 1:numel(newData)
    % Append the data:
    lineData = [lineData, newData(k).xLine, newData(k).yLine, ...
        newData(k).zLine, styleData];
    pointData = [pointData, newData(k).xPoints, newData(k).yPoints, ...
        newData(k).zPoints, styleData];
    jumpData = [jumpData, newData(k).xJumps, newData(k).yJumps, ...
        newData(k).zJumps, styleData];
end

%% Plotting starts here:

% Plot the curve
h1 = plot3(lineData{:});
set(h1, 'Marker', 'none', lineStyle{:})

% Ensure the plot is held:
hold on

% Plot the points:
h2 = plot3(pointData{:});
% Change the style accordingly:
set(h2, 'LineStyle', 'none', pointStyle{:})

% Plot the jumps:;
if ( isempty(jumpData) || ischar(jumpData{1}) )
    jumpData = {[]};
end
h3 = plot3(jumpData{:});
% Change the style accordingly:
if ( isempty(jumpStyle) )
    set(h3, 'LineStyle', ':', 'Marker', 'none')
else
    set(h3, jumpStyle{:});
end

% Plot the dummy data, which includes both line and point style:
hDummy = plot3(lineData{:}); % See CHEBFUN/PLOT() for details.
if ( ~isempty(lineStyle) || ~isempty(pointStyle) )
    set(hDummy, lineStyle{:}, pointStyle{:});
end

%% Misc:

% Return hold state to what it was before:
if ( ~holdState )
    hold off
end

% We don't want these guys to be included in LEGEND(), so we turn them off.
set(h1, 'handlevis', 'off');
set(h2, 'handlevis', 'off');
set(h3, 'handlevis', 'off');
% The dummy plot is invisible, but its handle is visible (for LEGEND).
set(hDummy, 'handlevis', 'on', 'visible', 'off');     %  ¯\_(o.O)_/¯

% Give an output to the plot handles if requested:
if ( nargout > 0 )
    varargout = {h1 ; h2 ; h3 ; hDummy};
end

%%

end

function f = expand(f, numCols)
if ( numel(f) ~= numCols )
    if ( numel(f) == 1 )
        f = repmat(f, 1, numCols);
    else
        error('CHEBFUN:CHEBFUN:plot3:expand:dim', ...
              'CHEBFUN objects must have the same number of columns.');
    end
end

end

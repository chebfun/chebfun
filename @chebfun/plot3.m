function varargout = plot3(f, g, h, varargin)
%PLOT3   Plot for CHEBFUN objects in 3-D space.
%   PLOT3() is a three-dimensional analogue of PLOT().
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

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% [TODO]: Implement plotting of delta functions.

% Deal with an empty input:
if ( isempty(f) || isempty(g) || isempty(h) )
    if ( nargout == 1 )
        varargout{1} = plot3([]);
    end
    return
end

% Store the hold state of the current axis:
holdState = ishold;

% We can only plot real against real:
if ( ~isreal(f) || ~isreal(g) || ~isreal(h) )
    warning('CHEBFUN:plot:complex', ...
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


% Deal with 'jumpLine' input.
[jumpStyle, varargin] = chebfun.parseJumpStyle(varargin{:});

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

% Plot the curve
h1 = plot3(lineData{:});
set(h1, 'Marker', 'none')

% Ensure the plot is held:
hold on

% Plot the points:
h2 = plot3(pointData{:});
% Change the style accordingly:
set(h2, 'LineStyle', 'none')

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

function f = expand(f, numCols)
if ( numel(f) ~= numCols )
    if ( numel(f) == 1 )
        f = repmat(f, 1, numCols);
    else
        error('CHEBFUN:plot:dim', ...
              'CHEBFUN objects must have the same number of columns.');
    end
end

end

function varargout = waterfall(f, varargin)
%WATERFALL Waterfall plot for array-valued CHEBFUN objects and quasimatrices.
%   WATERFALL(U), or WATERFALL(U, T) where LENGTH(T) = MIN(SIZE(U)), plots a
%   "waterfall" plot of an array-valued CHEBFUN or quasimatrix. Unlike the
%   standard Matlab waterfall, chebfun/waterfall does not fill in the column
%   planes with opaque whitespace or connect edges to zero.
%
%   Additional plotting options can also be passed, for example
%       WATERFALL(U,T,'linewidth',2).
%
% See also PLOT, PLOT3, PLOTDATA.

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

% Check for 2nd input:
if ( isnumeric(varargin{1}) )
    t = varargin{1};
    varargin(1) = [];
else
    t = 1:numColumns(f);
end

if ( ~isreal(f) || ~isreal(t) )
     warning('CHEBFUN:waterfall:complex', ...
                'Imaginary parts of complex X and/or Y arguments ignored.');
    f = real(f);
    t = real(t);
end

% Convert array-valued CHEBFUN to a quasimatrix:
f = cheb2quasi(f);

% Style data:
styleData = varargin;

% Loop over the columns:
lineData = {}; pointData = {}; jumpData = {};
for k = 1:numel(f)
    % Get the data from PLOTDATA():
    newData(k) = plotData(f(k));

    % Get the y data:
    ydl = repmat(t(k), size(newData(k).xLine, 1), 1);
    ydp = repmat(t(k), size(newData(k).xPoints, 1), 1);
    ydj = repmat(t(k), size(newData(k).xJumps, 1), 1);

    % Append new data:
    lineData = [lineData, newData(k).xLine, ydl, newData(k).yLine, styleData]; 
    pointData = [pointData, newData(k).xPoints, ydp, newData(k).yPoints, styleData];
    jumpData = [jumpData, newData(k).xJumps, ydj, newData(k).yJumps, styleData];
end
    
% Plot the lines:
h1 = plot3(lineData{:});
set(h1, 'Marker', 'none')

% Ensure the plot is held:
hold on

% Plot the points:
h2 = plot3(pointData{:});
% Change the style accordingly:
set(h2, 'LineStyle', 'none')

% Plot the jumps:
if ( isempty(jumpData) || ischar(jumpData{1}) )
    jumpData = {[]};
end
h3 = plot3(jumpData{:});
% Change the style accordingly:
set(h3, 'LineStyle', ':', 'Marker', 'none')

% Return hold state to what it was before:
if ( ~holdState )
    hold off
end

% Give an output to the plot handles if requested:
if ( nargout > 0 )
    varargout = {h1 ; h2 ; h3};
end

end
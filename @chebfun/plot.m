function varargout = plot(varargin)
%PLOT   Basic linear plot for CHEBFUN objects.
%   PLOT(F) plots the CHEBFUN object F in the interval where it is defined. If F
%   is complex valued, PLOT(F) is equivalent to PLOT(real(F), imag(F)).
%
%   PLOT(F, G) plots the CHEBFUN G versus the CHEBFUN F. Quasimatrices and
%   array-valued CHEBUFUN objects are also supported in the natural way.
%
%   PLOT(F, S) allows various line types, plot symbols, and colors to be used
%   when S is a character string made from one element from any or all the
%   following 3 columns:
%
%            b     blue          .     point              -     solid
%            g     green         o     circle             :     dotted
%            r     red           x     x-mark             -.    dashdot
%            c     cyan          +     plus               --    dashed
%            m     magenta       *     star             (none)  no line
%            y     yellow        s     square
%            k     black         d     diamond
%            w     white         v     triangle (down)
%                                ^     triangle (up)
%                                <     triangle (left)
%                                >     triangle (right)
%                                p     pentagram
%                                h     hexagram
%
%   The entries from the centre columns are plotted at the grid being used to
%   represent F (typically Chebyshev).
%
%   The X,Y pairs, or X,Y,S triples, can be followed by parameter/value pairs to
%   specify additional properties of the lines. For example,
%            f = chebfun(@sin);
%            plot(f, 'LineWidth', 2, 'Color', [.6 0 0])
%   will create a plot with a dark red line width of 2 points.
%
%   PLOT(F1, S1, F2, S2, F3, S3, ...) or PLOT(F1, G1, S1, F2, G2, S2, ...)
%   combines the plots defined by the (F,G,S) triples or (F,S) doubles, where
%   the F's and G's are CHEBFUN object and the S's are strings.
%
%   [HLINE, HPOINT, HJUMP] = PLOT(F) returns column vectors of handles to
%   lineseries objects (one for each column in the case of array-valued CHEBFUN
%   objects), corresponding to the line, point, and jump plots, respectively.
%
%   PLOT(F, 'interval', [A, B]) restricts the plot to the interval [A, B], which
%   can be useful when the domain of F is infinite, or for 'zooming in' on, say,
%   oscillatory CHEBFUN objects. If plotting an array-valued CHEBFUN or more
%   than one CHEBFUN in a call like PLOT(F, 'b', G, '--r', 'interval', [A, B])
%   this property is applied globally. Markers, such as 'o', or '.', are ignored
%   if the interval flag is used.
%
%   Besides the usual parameters that control the specifications of lines (see
%   linespec), the parameter JumpLine and DeltaLine determines the linestyle for
%   the discontinuities and the delta functions of the CHEBFUN F, respectively.
%   For example, PLOT(F, 'JumpLine', '-r') will plot discontinuities as solid
%   red lines and PLOT(F, 'deltaline', '-r') will plot the delta functions as
%   solid red lines. By default the plotting style for JumpLine is ':', and '-'
%   for delta functions and colours are chosen to match the lines they
%   correspond to. It is possible to modify other properties of JumpLines syntax
%   like PLOT(F, 'JumpLine', {'color', 'r', 'LineWidth', 5}). JumpLines and
%   deltaLines can be suppressed with the argument 'none'.
%
%   Note that the PLOT(F, 'numpts', N) option for V4 is deprecated, and this
%   call now has no effect.
%
% See also PLOTDATA, PLOT3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%  We actually plot a fifth 'dummy' plot which contains both the line and marker
%  styles. This is the only plot which has 'handlevis', 'on', and so will be
%  included in subsequent calls to LEGEND(). The handle to this plot is included
%  as the fifth output, and the plot itself is set to 'visible', 'off'. Note
%  that if a user modifies the plot by adjusting HLINE or HPOINT, subsequent
%  calls the LEGEND will not have the updated style. There's no way to get
%  around this without a complicated callback procedure.
%
%  Note, to clarify, the reason we need all this is that we must plot the line
%  and points as different plots. MATLAB/PLOT() doesn't have this problem.
%
%  Also note that the handles are now output as individually, rather than
%  forced in to an array.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: Document plotting of DELTFAFUN objects. (In particular, handle output)

%% Deal with an empty input:
if ( isempty(varargin{1}) )
    if ( nargout == 1 )
        varargout{1} = plot([]);
    end
    return
end

%% Initialization:

% Suppress inevitable warning for growing these arrays:
%#ok<*AGROW>

% Store the hold state of the current axis:
holdState = ishold;

% Store the current X and Y-limits, and whether ylim is in manual or auto mode.
if ( holdState )
    xLimCurrent = get(gca, 'xlim');
    yLimCurrent = get(gca, 'ylim');
    yLimModeCurrent = get(gca, 'ylimmode');
end

% Initialize flags:
isComplex = false;
intervalIsSet = false;

% Initialize XLIM and YLIM. Note that the first entries are initialized to be
% UPPER limits on the LOWER parts of the axes, while the second entries
% correspond to LOWER bounds on the UPPER parts of the axes. Hence, this
% somewhat strange convention.
xLim = [inf, -inf];
yLim = [inf, -inf];
defaultXLim = 1;
defaultYLim = 1;

% Check to see if the 'interval' flag has been set:
interval = [];
loc = find(strcmpi(varargin, 'interval'));
if ( any(loc) )
    interval = varargin{loc+1};
    varargin(loc:loc+1) = [];
    intervalIsSet = true;
else
    % TODO: Do we want to support this?
    % Allow plot(f, [a, b]) as shorthand for plot(f, 'interval', [a, b]):
    loc = find(cellfun(@(f) isa(f, 'chebfun'), varargin));
    for k = 1:numel(loc)
        if ( loc(k) < nargin && isnumeric(varargin{loc(k)+1}) )
            interval = varargin{loc(k)+1};
            varargin(loc(k)+1) = [];
            intervalIsSet = true;
            break
        end
    end
end

% Support 'interval' by evaluating on a fixed grid size (2000 points). See #602.
if ( intervalIsSet )
    for k = 1:numel(varargin)
        if ( isa(varargin{k}, 'chebfun') )
            dom = union(domain(varargin{k}), interval);
            dom(dom < interval(1) | dom > interval(end)) = [];
            varargin{k} = chebfun(@(x) feval(varargin{k}, x), dom, 2000);
        end
    end
end

% Initialise storage:
lineData = {};
pointData = {};
jumpData = {};
deltaData = {};
jumpLineIsSet = any(cellfun(@(v) ischar(v) && strcmpi(v, 'JumpLine'), varargin));

% Remove global plotting options from input arguments.
[lineStyle, pointStyle, jumpStyle, deltaStyle, varargin] = ...
    chebfun.parsePlotStyle(varargin{:});

%% Preparation of the data:

% Get the data for plotting from PLOTDATA():
while ( ~isempty(varargin) )

    % Acquire plotting data for each CHEBFUN / pair of CHEBFUNs:
    if ( (numel(varargin) > 1) && ... % PLOT(f, g).
            (isa(varargin{2}, 'chebfun') || isnumeric(varargin{2})) ) 
        % Remove CHEBFUN objects from array input:
        f = varargin{1};
        g = varargin{2};
        varargin(1:2) = [];
        
        % We can only plot real against real:
        isComplex = false;
        if ( ~isreal(f) || ~isreal(g) )
            warning('CHEBFUN:CHEBFUN:plot:complex', ...
                'Imaginary parts of complex X and/or Y arguments ignored.');
            f = real(f);
            g = real(g);
        end
        
        % Call PLOTDATA():
        if ( isnumeric(f) )
            % For plotting numeric data in a CHEBFUN/PLOT() call.
            newData = plotData(chebfun());
            newData.xLine = f;
            newData.yLine = g;
            newData.xPoints = f;
            newData.yPoints = g;
            newData.xJumps = NaN;
            newData.yJumps = NaN;  
            newData.xDeltas = NaN;
            % Do nothing
        elseif ( numel(f) == 1 && numel(g) == 1 )
            % Array-valued CHEBFUN case:
            newData = plotData(f, g);
        else
            % QUASIMATRIX case:
            f = num2cell(f);
            g = num2cell(g);
            if ( numel(f) > 1 && numel(g) > 1 )
                if ( numel(f) ~= numel(g) )
                    error('CHEBFUN:CHEBFUN:plot:dim', ...
                    'CHEBFUN objects must have the same number of columns.');
                end
                for k = 1:numel(f)
                    newData(k) = plotData(f{k}, g{k});
                end
            elseif ( numel(f) > 1 && numel(g) == 1 )
                for k = 1:numel(f)
                    newData(k) = plotData(f{k}, g{1});
                end
            elseif ( numel(f) == 1 && numel(g) > 1 )
                for k = 1:numel(f)
                    newData(k) = plotData(f{1}, g{k});
                end            
            end
        end

    else  % PLOT(f).
        
        % Remove CHEBFUN from array input:
        f = varargin{1};
        varargin(1) = [];
        isComplex = ~isreal(f);
        
        % Loop over the columns:
        for k = 1:numel(f)
            if ( isComplex )
                newData(k) = plotData(real(f(k)), imag(f(k)));
            else
                newData(k) = plotData(f(k));
            end
        end

    end
    
    % Style data.
    pos = 0;
    styleData = [];

    lv = length(varargin);
    % Find the location of the next CHEBFUN in the input array:
    while ( (pos < lv) && ~isa(varargin{pos+1}, 'chebfun') )
        if ( pos+1 < lv && isnumeric(varargin{pos+1}) && ...
                isnumeric(varargin{pos+2}) )
            break
        end
        pos = pos + 1;
    end
    
    if ( pos > 0 )
        styleData = varargin(1:pos);
        varargin(1:pos) = [];
        % Remove deprecated 'numpts' option:
        idx = find(strcmp(styleData, 'numpts'), 1);
        if ( any(idx) )
            styleData(idx:(idx+1)) = [];
        end
    end
    
    % Loop over the columns:
    for k = 1:numel(newData)
        
        % Update axis limits:
        xLim = [min(newData(k).xLim(1), xLim(1)), ...
            max(newData(k).xLim(2), xLim(2))];
        yLim = [min(newData(k).yLim(1), yLim(1)), ...
            max(newData(k).yLim(2), yLim(2))];

        % Append new data:
        lineData = [lineData, newData(k).xLine, newData(k).yLine, styleData];
        pointData = [pointData, newData(k).xPoints, newData(k).yPoints, ...
            styleData];
        jumpData = [jumpData, newData(k).xJumps, newData(k).yJumps, styleData];
        deltaData{k} = {newData(k).xDeltas, newData(k).yDeltas, ...
            newData(k).yDeltaBase, styleData};
        
        defaultXLim = defaultXLim & newData(k).defaultXLim;
        defaultYLim = defaultYLim & newData(k).defaultYLim;
    end
    
    % If xLim(1) == xLim(2), let Matlab figure out a proper xLim:
    if ( ~diff(xLim) )
        defaultXLim = 1;
    end
    
    % If yLim(1) == yLim(2), let Matlab figure out a proper yLim:
    if ( ~diff(yLim) )
        defaultYLim = 1;
    end
    
end

%% Plotting starts here:

% Acquire initial color cycle if running R2014b+.
if ( ~verLessThan('matlab', '8.4') )
    if ( ~holdState )
        set(gca, 'ColorOrderIndex', 1);
    end
    originalColorOrder = get(gca, 'ColorOrderIndex');
end

% Plot the lines:
h1 = plot(lineData{:});
set(h1, 'Marker', 'none', lineStyle{:})

% Ensure the plot is held:
hold on

% Get color cycle prior to point plot if running R2014b.
if ( ~verLessThan('matlab', '8.4') )
    newColorOrder = get(gca, 'ColorOrderIndex');
    set(gca, 'ColorOrderIndex', originalColorOrder)
end

% Plot the points:
h2 = plot(pointData{:});
% Change the style accordingly:
set(h2, 'LineStyle', 'none', pointStyle{:})
if ( intervalIsSet )
    % Markers are meaningless if the 'interval' flag is used.
    set(h2, 'Marker', 'none', pointStyle{:})
end

% Reset color cycle prior to jump plot if running R2014b.
if ( ~verLessThan('matlab', '8.4') )
    set(gca, 'ColorOrderIndex', originalColorOrder);
end

% Plot the jumps:
if ( isempty(jumpData) || ischar(jumpData{1}) )
    jumpData = {[]};
end
h3 = plot(jumpData{:});
% Change the style accordingly:
if ( ~isempty(jumpStyle) )
    set(h3, jumpStyle{:});
end
if ( ~jumpLineIsSet )
    if ( isComplex )
        set(h3, 'LineStyle', 'none') 
    else
        set(h3, 'LineStyle', ':')
    end
    set(h3, 'Marker', 'none') 
end

% Plot the Delta functions:
if ( isempty(deltaData) || ~isnumeric(deltaData{1}{1}) )
    h4 = plot([]);
else
    h4 = plotDeltas(deltaData);
end
if ( ~isempty(deltaStyle) )
    set(h4, deltaStyle{:});
end    

% Reset colors prior to legend data plot if running R2014b.
if ( ~verLessThan('matlab', '8.4') )
    set(gca, 'ColorOrderIndex', originalColorOrder);
end

% Generate data for the dummy plot for the legend.  dummyData is the same as
% lineData but with the numeric entries replaced by NaNs of suitable dimension
% so that they don't show up on the plot but are still "active" for the legend.
dummyData = lineData;
for i = 1:1:length(dummyData)
    if ( isnumeric(lineData{i}) )
        dummyData{i} = NaN(1, size(dummyData{i}, 2));
    end
end

% Plot the dummy data, which includes both line and point style:
hDummy = plot(dummyData{:});
if ( ~isempty(lineStyle) || ~isempty(pointStyle) )
    set(hDummy, lineStyle{:}, pointStyle{:});
end

% Reset colors prior to legend data plot if running R2014b.
if ( ~verLessThan('matlab', '8.4') )
    set(gca, 'ColorOrderIndex', newColorOrder);
end

%% Setting xLim and yLim:

% If holding, then make sure not to shrink the X-limits.
if ( holdState )
    xLim = [min(xLimCurrent(1), xLim(1)), max(xLimCurrent(2), xLim(2))];
    yLim = [min(yLimCurrent(1), yLim(1)), max(yLimCurrent(2), yLim(2))];
end

% We always want to set the x-limits. Otherwise, plots like
%   plot(chebfun(@(x) sin(x), [0 pi])
% will have extra white space around the ends of the domain, and look ugly.

if ( diff(xLim) )
    set(gca, 'xLim', xLim)
end

% Set the Y-limits if appropriate values have been suggested, or if we were
% holding on when we entered this method:
if ( ~defaultYLim || (holdState && strcmp(yLimModeCurrent, 'manual')) ) && ...
        ( diff(yLim) )
    set(gca, 'yLim', yLim)
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
set(h4, 'handlevis', 'off');

% Ensure that the both the handle for the dummy plot and the dummy plot itself
% are visible for LEGEND().  The dummy plot could in theory have its 'Visible'
% property set to 'off', but the data is NaN, so it doesn't draw anything
% anyway.  Moreover, if 'Visible' is 'off', the entries in the legend will be
% grayed out in R2016b and later.
set(hDummy, 'handlevis', 'on', 'visible', 'on');     %  ¯\_(o.O)_/¯

% Give an output to the plot handles if requested:
if ( nargout > 0 )
    varargout = {h1 ; h2 ; h3 ; h4 ; hDummy};
end

end

function h = plotDeltas(deltaData)
%PLOTDELTAS   Plots delta functions.
    h = [];

    % Get and save the current ColorOrder if running on R2014a or earlier.
    if ( verLessThan('matlab', '8.4') )
        originalColorOrder = get(gca, 'ColorOrder');
        colorOrder = circshift(originalColorOrder, 1);
    end

    for k = 1:1:numel(deltaData)
        % Set color for the next delta function plot.
        if ( verLessThan('matlab', '8.4') )
            % Manually manipulate the ColorOrder for R2014a or earlier.
            colorOrder = circshift(colorOrder, -1);
            set(gca, 'ColorOrder', colorOrder);
        else
            % Use ColorOrderIndex for R2014b and later.
            set(gca, 'ColorOrderIndex', k);
        end

        h = [h, plotDeltasForOneFunction(deltaData{k}{:})];
    end

    % Restore the ColorOrder if running on R2014a or earlier.
    if ( verLessThan('matlab', '8.4') )
        set(gca, 'ColorOrder', originalColorOrder);
    end
end

function h = plotDeltasForOneFunction(varargin)
%PLOTDELTASFORONEFUNCTION   Plots delta functions for one function.

h = [];
j = 1;
while ( ~isempty(varargin) )
    % Extract data which is in triplets:
    xData = varargin{1};  % locations of delta functions
    yData = varargin{2};  % magnitude of delta functions
    yBase = varargin{3};  % starting height for delta functions
    
    % Delete these arguments, since they have been copied:
    varargin(1:3) = [];

    % Handle the delta style argument:
    style = '';
    if ( ~isempty(varargin{1}) && iscell(varargin{1}) )
        style = varargin{1}{1};
    end
    varargin(1) = [];
    
    % Ignore complete NaN data:
    if ( all(isnan(xData)) )
        continue
    end
    
    if ( isnumeric(xData) )
        % Remove mixed NaN data:
        nanIdx = isnan(xData);
        xData(nanIdx) = [];
        yData(nanIdx) = [];
        yBase(nanIdx) = [];

        % merge duplicate delta functions.
        [yData, xData, dupIdx] = deltafun.mergeColumns(yData.', xData.');
        xData = xData.';
        yData = yData.';
        
        % For delta functions at discontinuities, set the base value to the
        % average of the function values:
        yBase(dupIdx-1) = 1/2*(yBase(dupIdx-1)+yBase(dupIdx));
        yBase(dupIdx) = [];
        yFinish = yBase + yData;
                        
        % Add NaNs in order to plot isolated vertical lines:

        % NaNs for xData:
        xFullData = zeros(3*size(xData,1), 1);
        xFullData(1:3:end) = xData;
        xFullData(2:3:end) = xData;
        xFullData(3:3:end) = NaN;
        
        % NaNs for yData
        yFullData = zeros(3*size(xData,1), 1);
        yFullData(1:3:end) = yBase;
        yFullData(2:3:end) = yFinish;
        yFullData(3:3:end) = NaN;
        
        fullData{1} = xFullData;
        fullData{2} = yFullData;        
        fullData{3} = style;
    end

    % Plot the vertical lines:
    h(j) = plot(fullData{:});
    
    % Plot the markers for delta functions:
    for jj = 1:length(xData)
        % For positive delta functions, we plot an arrow:
        if ( yData(jj) >= 0 )
            marker = '^';
        end
        % For negative delta functions, we plot a v: 
        if ( yData(jj) < 0 )
            marker = 'v';
        end        

        % Plot the markers and match the color with the delta lines:
        hjj = plot(xData(jj), yFinish(jj), marker, ...
            'markersize', 6, 'linestyle', 'none', 'handlevis', 'off');
        set(hjj, 'color', get(h(j), 'color'));
        set(hjj, 'markerfacecolor', get(h(j), 'color'))
    end
        
    j = j + 1;
end

end

function varargout = plot(varargin)
%PLOT   Basic linear plot for CHEBFUN objects.
%   PLOT(F) plots the CHEBFUN object F.
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
%   this property is applied globally.
%
%   Besides the usual parameters that control the specifications of lines (see
%   linespec), the parameter JumpLines determines the linestyle for
%   discontinuities of the CHEBFUN F. For example, PLOT(F, 'JumpLine', '-r')
%   will plot discontinuities as solid red lines. By default the plotting style
%   is ':', and colours are chosen to match the lines they correspond to. It is
%   possible to modify other properties of JumpLines syntax like PLOT(F,
%   'JumpLine', {'color', 'r', 'LineWidth', 5}). JumpLines can be suppressed
%   with the argument 'JumpLine','none'.
%
%   Note that the PLOT(F, 'numpts', N) option for V4 is deprecated, and this
%   call now has no effect.
%
% See also PLOTDATA, PLOT3.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% TODO: Figure out the y axis limit for functions which blow up.

% Deal with an empty input:
if ( isempty(varargin{1}) )
    if ( nargout == 1 )
        varargout{1} = plot([]);
    end
    return
end

% Store the hold state of the current axis:
holdState = ishold;

% Store the current X and Y-limits:
if ( holdState )
    xLimCurrent = get(gca, 'xlim');
    yLimCurrent = get(gca, 'ylim');
end

% Initialize flags:
isComplex = false;
intervalIsSet = false;
xLim = [inf, -inf];
yLim = [inf, -inf];

% Suppress inevitable warning for growing these arrays:
%#ok<*AGROW>

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

% Initialise storage:
lineData = {};
pointData = {};
jumpData = {};
deltaData = {};

% Remove global plotting options from input arguments.
[lineStyle, pointStyle, jumpStyle, varargin] = ...
    chebfun.parsePlotStyle(varargin{:});

%%
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
            warning('CHEBFUN:plot:complex', ...
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
            newData.xDeltas = [];
            newData.yDeltas = [];
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
                    error('CHEBFUN:plot:dim', ...
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

    else                                                       % PLOT(f).
        
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
        
        % Handle the 'interval' flag:
        if ( ~isComplex && intervalIsSet && (size(newData(k).xLine, 2) == 1) )
            ind = newData(k).xLine < interval(1) | ...
                newData(k).xLine > interval(end);
            newData(k).xLine(ind) = [];
            newData(k).yLine(ind,:) = [];
            ind = newData(k).xPoints < interval(1) | ...
                newData(k).xPoints > interval(end);
            newData(k).xPoints(ind) = [];
            newData(k).yPoints(ind,:) = [];
            ind = newData(k).xJumps < interval(1) | ...
                newData(k).xJumps > interval(end);
            newData(k).xJumps(ind) = [];
            newData(k).yJumps(ind,:) = [];            
            ind = newData(k).xDeltas < interval(1) | ...
                newData(k).xDeltas > interval(end);
            newData(k).xDeltas(ind) = [];
            newData(k).yDeltas(ind,:) = [];
            
            newData(k).xLim = interval;            
        end
        
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
        deltaData = [deltaData, newData(k).xDeltas, newData(k).yDeltas];
    end
    
    % If xLim(1) == xLim(2), set xLim [inf -inf] and let Matlab figure out a
    % proper xLim:
    if ( ~diff(xLim) )
        xLim = [inf, -inf];
    end
    
    % If yLim(1) == yLim(2), set yLim [inf -inf] and let Matlab figure out a
    % proper yLim:
    if ( ~diff(yLim) )
        yLim = [inf, -inf];
    end
    
end
% Plot the lines:
h1 = plot(lineData{:});
set(h1, 'Marker', 'none', lineStyle{:})

% Ensure the plot is held:
hold on

% Plot the points:
h2 = plot(pointData{:});
% Change the style accordingly:
set(h2, 'LineStyle', 'none', pointStyle{:})

% Plot the jumps:
if ( isempty(jumpData) || ischar(jumpData{1}) )
    jumpData = {[]};
end
h3 = plot(jumpData{:});
% Change the style accordingly:
if ( isempty(jumpStyle) )
    if ( isComplex )
        %[TODO]: The following statement can not be reached:
        set(h3, 'LineStyle', 'none', 'Marker', 'none')
    else
        set(h3, 'LineStyle', ':', 'Marker', 'none')
    end
else
    set(h3, jumpStyle{:});
end

% Plot the Delta functions:
if ( isempty(deltaData) )
    deltaData = {[]};
end
h4 = stem(deltaData{:}, 'd', 'fill');

%% 
% Do we want a style for delta functions?
% if ( isempty(jumpStyle) )
%     if ( isComplex )
%         %[TODO]: The following statement can not be reached:
%         set(h3, 'LineStyle', 'none', 'Marker', 'none')
%     else
%         set(h3, 'LineStyle', ':', 'Marker', 'none')
%     end
% else
%     set(h3, jumpStyle{:});
% end
% Set the X-limits if appropriate values have been suggested:
if ( all(isfinite(xLim)) )

    % If holding, then make sure not to shrink the X-limits.
    if ( holdState )
        xLim = [min(xLimCurrent(1), xLim(1)), max(xLimCurrent(2), xLim(2))];
    end

    set(gca, 'xlim', sort(xLim))
end

% Set the Y-limits if appropriate values have been suggested:
if ( all(isfinite(yLim)) )

    % If holding, then make sure not to shrink the Y-limits.
    if ( holdState )
        yLim = [min(yLimCurrent(1), yLim(1)), max(yLimCurrent(2), yLim(2))];
    end

    set(gca, 'ylim', sort(yLim))
end

% Return hold state to what it was before:
if ( ~holdState )
    hold off
end

% Give an output to the plot handles if requested:
if ( nargout > 0 )
    varargout = {h1 ; h2 ; h3 ; h4};
end

end


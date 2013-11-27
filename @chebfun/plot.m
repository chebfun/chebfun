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
% See also PLOTDATA, PLOT3.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% [TODO]: Implement plotting of delta functions.

% Deal with an empty input:
if ( isempty(varargin{1}) )
    if ( nargout == 1 )
        varargout{1} = plot([]);
    end
    return
end

% Store the hold state of the current axis:
holdState = ishold;

% Initialise storage:
lineData = {};
pointData = {};
jumpData = {};
isComplex = false;

% Suppress inevitable warning for growing these arrays:
%#ok<*AGROW>

% Check to see if the 'interval' flag has been set:
intervalIsSet = 0;
loc = find(strcmpi(varargin, 'interval'));
if ( any(loc) )
    intervalIsSet = 1;
    interval = varargin{loc+1};
    varargin(loc:loc+1) = [];
else
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
    
%%
% Get the data for plotting from PLOTDATA():
while ( ~isempty(varargin) ) 

    % Acquire plotting data for each CHEBFUN / pair of CHEBFUNs:
    if ( (numel(varargin) > 1) && isa(varargin{2}, 'chebfun') ) % PLOT(f, g).
        f = varargin{1};
        g = varargin{2};
        
        % We can only plot real against real:
        if ( ~isreal(f) || ~isreal(g) )
            warning('CHEBFUN:plot:complex', ...
                'Imaginary parts of complex X and/or Y arguments ignored.');
            f = real(f);
            g = real(g);
        end
        
        % Call PLOTDATA():
        newData = plotData(f, g);
        % Remove CHEBFUN objects from array input:
        varargin(1:2) = [];
        
    else                                                       % PLOT(f).
        % Call PLOTDATA():
        newData = plotData(varargin{1});
        % Remove CHEBFUN from array input:
        varargin(1) = [];
        
    end
    
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
    
    if ( ~isreal( newData.yLine ) ) % Deal with complex-valued functions.
        % Assign x to be the real part, and y to be the imagiary part:
        newData.xLine = real(newData.yLine);
        newData.yLine = imag(newData.yLine);
        newData.xPoints = real(newData.yPoints);
        newData.yPoints = imag(newData.yPoints);
        newData.xJumps = real(newData.yJumps);
        newData.yJumps = imag(newData.yJumps);
        isComplex = true;
    elseif ( intervalIsSet && size(newData.xLine, 2) == 1 ) % Deal with 'interval' flag.
        idx = newData.xLine < interval(1) | newData.xLine > interval(end);
        newData.xLine(idx) = [];
        newData.yLine(idx,:) = [];
        idx = newData.xPoints < interval(1) | newData.xPoints > interval(end);
        newData.xPoints(idx) = [];
        newData.yPoints(idx,:) = [];
        idx = newData.xJumps < interval(1) | newData.xJumps > interval(end);
        newData.xJumps(idx) = [];
        newData.yJumps(idx,:) = [];
    end
    
    % Append new data to the arrays which will be passed to built in PLOT():
    lineData = [lineData, newData.xLine, newData.yLine, styleData]; 
    pointData = [pointData, newData.xPoints, newData.yPoints, styleData];
    jumpData = [jumpData, newData.xJumps, newData.yJumps, styleData];
     
end

% Plot the lines:
h1 = plot(lineData{:});
set(h1, 'Marker', 'none')

% Ensure the plot is held:
hold on

% Plot the points:
h2 = plot(pointData{:});
% Change the style accordingly:
set(h2, 'LineStyle', 'none')

% Plot the jumps:
if ( isempty(jumpData) )
    jumpData = {[]};
end
h3 = plot(jumpData{:});
% Change the style accordingly:
if ( isComplex )
    set(h3, 'LineStyle', 'none', 'Marker', 'none')
else
    set(h3, 'LineStyle', ':', 'Marker', 'none')
end

% Return hold state to what it was before:
if ( ~holdState )
    hold off
end

% Give an output to the plot handles if requested:
if ( nargout > 0 )
    varargout = {h1 ; h2 ; h3};
end

end

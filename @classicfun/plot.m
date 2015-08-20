function varargout = plot(f, varargin)
%PLOT   Basic linear plot for CLASSICFUN objects. 
%   PLOT(F) plots the CLASSICFUN object F.
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
%   The X,Y pairs, or X,Y,S triples, can be followed by parameter/value
%   pairs to specify additional properties of the lines. For example,
%            F = bndfun(@sin, [-2 7]);
%            plot(F, 'LineWidth', 2, 'Color', [.6 0 0]) 
%   will create a plot with a dark red line of width 2 points.
%
%   The line data for the plot is determined by PLOTDATA, which queries F.onefun
%   for suitable information. CLASSICFUN/PLOT will display points at sample values (if
%   returned by PLOTDATA) on if LENGTH(F) < 257.
%
%   H1 = PLOT(F) returns a column vector of handles to lineseries objects, one
%   handle per plotted line (in the case of array-valued CLASSICFUN objects).    
%   [H1, H2] returns a second vector of column handles, this time for each of
%   the marker plots.
%
% See also PLOTDATA.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Deal with an empty input:
if ( isempty(f) )
    if ( nargout == 1 )
        varargout{1} = plot([]);
    end
    return
end

% Store the hold state of the current axis:
holdState = ishold;

%%
% Get the data for plotting from PLOTDATA():
if ( nargin > 1 && isa(varargin{1}, 'classicfun') )
    % Deal with plot(f, g);
    g = varargin{1};
    varargin(1) = [];
    % We can only plot real against real:
    f = real(f); 
    g = real(g);
    % Call PLOTDATA():
    data = plotData(f, g);
else
    g = [];
    data = plotData(f);
end

%%
% Classicfun has no delta functions:
data.xDeltas = [];
data.yDeltas = [];

%%
% Plot the curve:
if ( isreal(f) )
    h1 = plot(data.xLine, data.yLine, varargin{:}); 
else
    h1 = plot(data.yLine, varargin{:}); 
end

% No points on the line plot:
set(h1, 'Marker', 'none') 

% Hold the plot:
hold on

%%
% Plot the points:
if ( isreal(f) )
    h2 = plot(data.xPoints, data.yPoints, varargin{:});
else
    h2 = plot(data.yPoints, varargin{:});
end

% Change the style accordingly:
set(h2, 'LineStyle', 'none')
if ( all(strcmp(get(h2, 'Marker'), 'none')) && length(f) < 257 )
    set(h2,'Marker', 'o')
end

set(gca, 'xlim', data.xLim)
set(gca, 'ylim', data.yLim)

%%
% Return hold state to what it was before:
if ( ~holdState )
    hold off
end

% Give an output if one was requested:
if ( nargout > 0 )
    varargout{1} = h1;
    varargout{2} = h2;
end

end

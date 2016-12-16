function varargout = plot(f, varargin)
%PLOT   Basic linear plot for TRIGTECH objects. 
%   PLOT(F) plots the TRIGTECH object F.
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
%   The entries from the centre columns are plotted at TRIGPTS being
%   used to represent F. If no options from this column are chosen, 'o' is
%   chosen by default if LENGTH(F) < 256.
%
%   The X,Y pairs, or X,Y,S triples, can be followed by parameter/value pairs to
%   specify additional properties of the lines. For example, 
%            f = trigtech.constructor(@sin);
%            plot(f, 'LineWidth', 2, 'Color', [.6 0 0]) 
%   will create a plot with a dark red line width of 2 points.
%
%   H1 = PLOT(F) returns a column vector of handles to lineseries objects, one
%   handle per plotted line (in the case of array-valued TRIGTECH objects).
%   [H1, H2] returns a second vector of column handles, this time for each of
%   the marker plots.
%
% See also PLOT3, PLOTDATA.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
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
if ( nargin > 1 && isa(varargin{1}, 'trigtech') )
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
% Plot the curve:
if ( isreal(f) )
    h1 = plot(real(data.xLine), real(data.yLine), varargin{:}); 
else
    h1 = plot(data.yLine, varargin{:}); 
end
set(h1, 'Marker', 'none') 
hold on

%%
% Plot the points:
if ( isreal(f) )
    h2 = plot(real(data.xPoints), real(data.yPoints), varargin{:});
else
    h2 = plot(data.yPoints, varargin{:});
end

% Change the style accordingly:
set(h2,'LineStyle', 'none')
if ( all(strcmp(get(h2, 'Marker'), 'none')) && (length(f) < 257) )
    set(h2, 'Marker', 'o')
end

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

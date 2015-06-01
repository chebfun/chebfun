function varargout = plot3(f, g, h, varargin)
%PLOT3   Plot for FUN objects in 3-D space. 
%   PLOT3() is a three-dimensional analogue of PLOT().
%   
%   PLOT3(X, Y, Z), where Z, Y, and Z are three CHEBTECH objects, plots a line
%   in 3-space. X, Y, and Z may be array-valued, but must have the same number
%   of columns.
%   
%   Various line types, plot symbols, and colors may be obtained with PLOT3(X,
%   Y, Z, S) where s is a 1, 2 or 3 character string made from the characters
%   listed under the PLOT command.
%
%   H1 = PLOT3(X, Y, Z) returns a column vector of handles to lineseries
%   objects, one handle per plotted line (in the case of vector-valued CHEBTECH
%   objects). [H1, H2] returns a second vector of column handles, this time for
%   each of the marker plots.
%
% See also PLOT, PLOTDATA.

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
% We can only plot real data in this way:
f = real(f);
g = real(g);
h = real(h);

% Get the data for plotting from PLOTDATA():
data = plotData(f, g, h);

%%
% Plot the curve:
h1 = plot3(data.xLine, data.yLine, data.zLine); 
set(h1, 'Marker', 'none') 
hold on

%%
% Plot the points:
h2 = plot3(data.xPoints, data.yPoints, data.zPoints, varargin{:});

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

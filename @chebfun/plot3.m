function varargout = plot3(f, g, h, varargin)
%PLOT3   Plot for CHEBFUN objects in 3-D space. 
%   PLOT() is a three-dimensional analogue of PLOT().
%   
%   PLOT3(X, Y, Z), where X, Y, and Z are three CHEBFUN objects, plots a line in
%   3-space. X, Y, and Z may be array-valued, but must have the same number of
%   columns.
%   
%   Various line types, plot symbols, and colors may be obtained with PLOT3(X,
%   Y, Z, S) where S is a string of length 1, 2 or 3 containing characters
%   listed under the PLOT command.
%
%   H1 = PLOT3(F, G, H) returns a column vector of handles to lineseries
%   objects, one handle per plotted line (in the case of vector-valued CHEBFUN
%   objects). [H1, H2] returns a second vector of column handles, this time for
%   each of the marker plots.
%  
%   Example: A helix:
%       
%       t = chebfun('t', [0, 10*pi]);
%       plot3(sin(t), cos(t), t);
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
% We can only plot real against real:
% [TODO]: Replace this once we have CHEBFUN/REAL() and CHEBFUN/IMAG().
% if ( ~isreal(f) || ~isreal(g) || ~isreal(h) )
%     warning('CHEBFUN:plot:complex', ...
%        'Warning: Imaginary parts of complex X and/or Y arguments ignored.');
% end
% f = real(f);
% g = real(g);
% h = real(h);

% Get the data for plotting from PLOTDATA():
data = plotData(f, g, h);

%%
% Plot the curve

h1 = plot3(data.xLine, data.yLine, data.zLine, varargin{:});
set(h1, 'Marker', 'none')
hold on

%%
% Plot the points:
h2 = plot3(data.xPoints, data.yPoints, data.zPoints, varargin{:});

% Change the style accordingly:
set(h2,'LineStyle', 'none')
if ( all(strcmp(get(h2, 'Marker'),'none')) ) && length(f) < 257
    set(h2,'Marker', 'none')
end

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
end

end
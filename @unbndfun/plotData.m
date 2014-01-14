function data = plotData(f, g)
%PLOTDATA    Useful data values for plotting a UNBNDFUN object.
%   DATA = PLOTDATA(F) returns a cell array of data values that can be used for
%   plotting F. In particular, DATA is a 4x1 cell array of the form {xLine,
%   fLine, xPoints, fPoints}, where xLine-fLine are a data pair for plotting the
%   continuous function F and xPoints-fPoints are the data pair for plotting
%   values of F on the underlying Chebyshev grid.
%
% See also PLOT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 || isempty(g) )
    % Get the data from the ONEFUN:
    data = plotData(f.onefun);
    
    % Map the 'x' data using f.mapping.for:
    data.xLine = f.mapping.for(data.xLine);
    data.xPoints = f.mapping.for(data.xPoints);
    
    % The plotting window is [-window, window] for doubly infinite domains, and
    % of width 2*window for singly infinite domains.
    window = 10;
    dom = f.domain;
    if ( all(isinf(dom)) ) % [-Inf, Inf]
        
        % Line data:
        data.yLine(data.xLine > window | data.xLine < -window) = [];
        data.xLine(data.xLine > window | data.xLine < -window) = [];
        % Point data:
        data.yPoints(data.xPoints > window | data.xPoints < -window) = [];
        data.xPoints(data.xPoints > window | data.xPoints < -window) = [];
        
    elseif ( isinf(dom(1)) ) % [-Inf, b]
        
        % Line data:
        data.yLine(data.xLine < dom(2) - 2*window) = [];
        data.xLine(data.xLine < dom(2) - 2*window) = [];
        % Point data:
        data.yPoints(data.xPoints < dom(2) - 2*window) = [];
        data.xPoints(data.xPoints < dom(2) - 2*window) = [];
        
    elseif ( isinf(dom(2)) ) % [a, Inf]
        
        % Line data:
        data.yLine(data.xLine > dom(1) + 2*window) = [];
        data.xLine(data.xLine > dom(1) + 2*window) = [];
        % Point data:
        data.yPoints(data.xPoints > dom(1) + 2*window) = [];
        data.xPoints(data.xPoints > dom(1) + 2*window) = [];
        
    end
    
end

end
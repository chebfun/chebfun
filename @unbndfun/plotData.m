function data = plotData(f)
%PLOTDATA    Useful data values for plotting a UNBNDFUN object.
%   DATA = PLOTDATA(F) returns a cell array of data values that can be used for
%   plotting F. In particular, DATA is a 4x1 cell array of the form {xLine,
%   fLine, xPoints, fPoints}, where xLine-fLine are a data pair for plotting the
%   continuous function F and xPoints-fPoints are the data pair for plotting
%   values of F on the underlying Chebyshev grid.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the data from the ONEFUN:
data = plotData(f.onefun);

% Map the 'x' data using f.mapping.for:
data{1} = f.mapping.for(data{1});
data{3} = f.mapping.for(data{3});

% The plotting window is [-window, window] for doubly infinite domains, and
% of width 2*window for singly infinite domains.
window = 10;
dom = f.domain;
if ( all(isinf(dom)) ) % [-Inf, Inf]
    
    % Line data:
    data{2}(data{1} > window | data{1} < -window) = [];
    data{1}(data{1} > window | data{1} < -window) = [];
    % Point data:
    data{4}(data{3} > window | data{3} < -window) = [];
    data{3}(data{3} > window | data{3} < -window) = [];
    
elseif ( isinf(dom(1)) ) % [-Inf, b]
    
    % Line data:
    data{2}(data{1} < dom(2) - 2*window) = [];
    data{1}(data{1} < dom(2) - 2*window) = [];
    % Point data:
    data{4}(data{3} < dom(2) - 2*window) = [];
    data{3}(data{3} < dom(2) - 2*window) = [];
    
elseif ( isinf(dom(2)) ) % [a, Inf]
    
    % Line data:
    data{2}(data{1} > dom(1) + 2*window) = [];
    data{1}(data{1} > dom(1) + 2*window) = [];
    % Point data:
    data{4}(data{3} > dom(1) + 2*window) = [];
    data{3}(data{3} > dom(1) + 2*window) = [];    
    
end
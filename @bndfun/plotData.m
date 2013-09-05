function data = plotData(f, g, h)
%PLOTDATA   Useful data values for plotting a BNDFUN object.
%   DATA = PLOTDATA(F) returns a struct containing data that can be used for
%   plotting F. In particular, DATA.xLine and DATA.yLine are for plotting smooth
%   curves (usually passed to plot with '-') and DATA.xPoints and DATA.yPoints
%   contain the (x, F(x)) data used to represent F.
%
%   DATA = PLOTDATA(F, G) returns data for PLOT(F, G), i.e., (F(x), G(x)), and
%   DATA = PLOTDATA(F, G, H) returns data for plots of the form PLOT3(F, G, H).
%   In the latter case, DATA also contains the fields zLine and zPoints.
%
% See also PLOT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the data from the ONEFUN:
if ( nargin == 1 || isempty(g) )
    % PLOT(F):
    data = plotData(f.onefun);
    % Map the 'x' data using f.mapping.for:
    data.xLine = f.mapping.for(data.xLine);
    data.xPoints = f.mapping.for(data.xPoints);
    
elseif ( nargin == 2 )
    % PLOT(F, G):
    data = plotData(f.onefun, g.onefun);
    
else
    % PLOT(F, G, H):
    data = plotData(f.onefun, g.onefun, h.onefun);
    
end

end
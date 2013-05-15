function data = plotData(f)
%PLOTDATA    Useful data values for plotting a BNDFUN object.
%   OUT = PLOTDATA(F) returns a struct containing data that can be used for
%   plotting F. In particular, DATA.xLine and DATA.fLine are for plotting smooth
%   curves (usually passed to plot with '-') and DATA.xPoints and DATA.yPoints
%   contain the (x, F(x)) data used to represent F.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the data from the ONEFUN:
data = plotData(f.onefun);

% Map the 'x' data using f.mapping.for:
data.xLine = f.mapping.for(data.xLine);
data.xPoints = f.mapping.for(data.xPoints);

end
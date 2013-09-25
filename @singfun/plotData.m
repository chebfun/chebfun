function data = plotData(f)
%PLOTDATA   Useful data values for plotting a SINGFUN object.
%   DATA = PLOTDATA(F) extracts PLOTDATA of the smooth part of F
%   and then scales it by the singular factors given in the EXPONENTS of F
%
% See also PLOT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%
% Get plot data from the smooth part:
data = plotData( f.smoothPart );
% Update extrapolated y-data:
x = data.xLine;
y = data.yLine;
y = y.*(x + 1).^f.exponents(1);
y = y.*(1 - x).^f.exponents(2);
data.yLine = y;

% Update sample point y-data:
y = data.yPoints;
x = data.xPoints;
y = y.*(x + 1).^f.exponents(1);
y = y.*(1 - x).^f.exponents(2);
data.yPoints = y;

end

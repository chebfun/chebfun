function data = plotData(f)
%PLOTDATA    Useful data values for plotting a SINGFUN object.
%   DATA = PLOTDATA(F) extracts PLOTDATA of the smooth part of F
%   and then scales it by the singular factors given in the EXPONENTS of F
%
% See also SMOOTHFUN/PLOTDATA.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%
% [TODO]: The following code may be simplified by calling SINGFUN.FEVAL()
% but what should the x-data be? (CHEBTECH uses 501 chebpts).
data = plotData( f.smoothPart );
% update extrapolated y-data 
x = data.xLine;
y = data.fLine;
y = y.*(x+1).^f.exponents(1);
y = y.*(1-x).^f.exponents(2);
data.fLine = y;

% update sample point y-data
y = data.fPoints;
x = data.xPoints;
y = y.*(x+1).^f.exponents(1);
y = y.*(1-x).^f.exponents(2);
data.fPoints = y;

end
function data = plotData(f)
%PLOTDATA    Useful data values for plotting a CHEBTECH object.
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

end
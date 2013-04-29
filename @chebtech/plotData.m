function data = plotData(f)
%PLOTDATA    Useful data values for plotting a CHEBTECH object.
%   DATA = PLOTDATA(F) returns a cell array of data values that can be used for
%   plotting F. In particular, DATA is a 4x1 cell array of the form {xLine,
%   fLine, xPoints, fPoints}, where xLine-fLine are a data pair for plotting the
%   continuous function F and xPoints-fPoints are the data pair for plotting
%   values of F on the underlying Chebyshev grid.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the number of points: (Oversample the wavelength)
npts = max(2001, round(4*pi*length(f)));

% Values on oversampled Chebyshev grid (faster than evaluating a uniform grid!).
xLine = f.chebpts(npts);
fLine = get(prolong(f, npts), 'values');

% Values on the Cheyshev grid tied to the CHEBTECH F:
xPoints = f.points();
fPoints = f.values;

% Wrap the data in a cell array:
data = {xLine, fLine, xPoints, fPoints};

end
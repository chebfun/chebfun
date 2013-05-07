function data = plotData(f, pref)
%PLOTDATA    Useful data values for plotting a CHEBTECH object.
%   OUT = PLOTDATA(F) returns a struct containing data that can be used for
%   plotting F. In particular, DATA.xLine and DATA.fLine are for plotting smooth
%   curves (usually passed to plot with '-') and DATA.xPoints and DATA.yPoints
%   contain the (x, F(x)) data at the Chebyshev grid F.POINTS() used to
%   represent F.
%
% See also PLOT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    pref = chebtech.pref();
end

% Get the number of points: (Oversample the wavelength)
npts = min(max(501, round(4*pi*length(f))), pref.chebtech.maxSamples);

% Initialise the output structure:
data = struct('xLine', [], 'fLine', [], 'xPoints', [], 'fPoints', []);

% Values on oversampled Chebyshev grid (faster than evaluating a uniform grid!).
data.xLine = f.chebpts(npts);
data.fLine = get(prolong(f, npts), 'values');

% Values on the Cheyshev grid tied to the CHEBTECH F:
data.xPoints = f.points();
data.fPoints = f.values;

end
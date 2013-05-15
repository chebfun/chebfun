function data = plotData(f, pref)
%PLOTDATA    Useful data values for plotting a CHEBTECH object.
%   DATA = PLOTDATA(F) returns a struct containing data that can be used for
%   plotting F. The struct DATA contains the following fields:
%
%       xLine: x-coordinates of for plotting smooth curves.
%       fLine: Function values of F at the coordinates stored in xLine.
%       xPoints: x-coordinates of the Chebyshev points used to represent F.
%       fPoints: Function value at the Chebyshev points used to represent F.
%   
%   DATA.xLine and DATA.fLine are used for plotting smooth curves (usually
%   passed to PLOT() with the '-' option). 
%
%   DATA.xPoints and DATA.fPoints contain the (x, F(x)) data at the Chebyshev
%   grid used to represent F, and are used for plots with markes (e.g.
%   PLOT(F,'-o').
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
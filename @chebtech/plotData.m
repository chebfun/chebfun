function data = plotData(f, g)
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
%   DATA = PLOTDATA(F, G) is similar but for plot calls of the form PLOT(F, G),
%   where both F and G are CHEBTECH objects.
%
% See also PLOT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    g = [];
end

% Get the number of points: (Oversample the wavelength)
len = max(length(f), length(g));
npts = min(max(501, round(4*pi*len)), chebtech.pref('maxSamples'));

% Initialise the output structure:
data = struct('xLine', [], 'fLine', [], 'xPoints', [], 'fPoints', []);
if ( isa(g, 'chebtech') )   % plot(f, g)
    % Values on oversampled Chebyshev grid (faster than evaluating a uniform grid!).
    data.xLine = get(prolong(f, npts), 'values');
    data.fLine = get(prolong(g, npts), 'values');

    % Values on the largest Cheyshev grid tied to the CHEBTECH objects F and G:
    data.xPoints = get(prolong(f, len), 'values');
    data.fPoints = get(prolong(g, len), 'values');
    
elseif ( isempty(g) )       % plot(f)
    % Values on oversampled Chebyshev grid (faster than evaluating a uniform grid!).
    data.xLine = f.chebpts(npts);
    data.fLine = get(prolong(f, npts), 'values');

    % Values on the Cheyshev grid tied to the CHEBTECH F:
    data.xPoints = f.points();
    data.fPoints = f.values;
    
else
    error('CHEBFUN:CHEBTECH:plotdata:DataType', ...
        'Cannot plot a %s against a %s.', class(f), calss(g));
    
end

end

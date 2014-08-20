function data = plotData(f, g, h)
%PLOTDATA   Useful data values for plotting a CHEBTECH object.
%   DATA = PLOTDATA(F) returns a struct containing data that can be used for
%   plotting F. The struct DATA contains the following fields:
%
%       xLine: x-coordinates of for plotting smooth curves.
%       yLine: Function values of F at the coordinates stored in xLine.
%       xPoints: x-coordinates of the Chebyshev points used to represent F.
%       yPoints: Function value at the Chebyshev points used to represent F.
%   
%   DATA.xLine and DATA.yLine are used for plotting smooth curves (usually
%   passed to PLOT() with the '-' option). 
%
%   DATA.xPoints and DATA.yPoints contain the (x, F(x)) data at the Chebyshev
%   grid used to represent F, and are used for plots with marks (e.g.
%   PLOT(F,'-o').
%
%   DATA = PLOTDATA(F, G) is similar but for plot calls of the form PLOT(F, G),
%   where both F and G are CHEBTECH objects. 
% 
%   DATA = PLOTDATA(F, G, H) is for plots of the form PLOT3(F, G, H). In this
%   instance, DATA also contains fields zLine and zPoints for the data
%   corresponding to H.
%
% See also PLOT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    g = [];
end
if ( nargin < 3 )
    h = [];
end

% Get the number of points: (Oversample the wavelength)
len = max([length(f), length(g), length(h)]);
npts = min(max(501, round(4*pi*len)), chebtech.techPref().maxLength);

% Initialise the output structure:
data = struct('xLine', [], 'yLine', [], 'xPoints', [], 'yPoints', [], ...
    'yLim', [], 'defaultXLim', 1, 'defaultYLim', 1);

if ( isempty(g) )       
    % PLOT(F):
    
    % Values on oversampled Chebyshev grid (faster than evaluating a uniform grid!).
    data.xLine = f.chebpts(npts);
    data.yLine = get(prolong(f, npts), 'values');

    % Values on the Cheyshev grid tied to the CHEBTECH F:
    data.xPoints = f.points();
    data.yPoints = f.coeffs2vals(f.coeffs);
    
    % yLim:
    data.yLim = [min(data.yLine(:)) max(data.yLine(:))];

elseif ( isa(g, 'chebtech') )   
    % PLOT(F, G)
    
    % Also return the grid points used.
    % Grid data for f:
    data.fGrid.xLine = f.chebpts(npts);
    % Use the maximum of the lenghts of f, g and h to match the number of
    % values returned:
    data.fGrid.xPoints = f.chebpts(len);
    
    % Grid data for g:
    data.gGrid.xLine = g.chebpts(npts);
    % Use the maximum of the lenghts of f, g and h to match the number of
    % values returned:    
    data.gGrid.xPoints = g.chebpts(len);
    
    % Values on oversampled Chebyshev grid (faster than evaluating a uniform grid!).
    data.xLine = get(prolong(f, npts), 'values');
    data.yLine = get(prolong(g, npts), 'values');

    % Values on the largest Cheyshev grid tied to the CHEBTECH objects F and G:
    data.xPoints = get(prolong(f, len), 'values');
    data.yPoints = get(prolong(g, len), 'values');
    
    % xLim:
    xdata = [get(f, 'lval'); data.xLine; get(f, 'rval')];
    data.xLim = [min(xdata(:)) max(xdata(:))];
    
    % yLim:
    ydata = [get(g, 'lval'); data.yLine; get(g, 'rval')];
    data.yLim = [min(ydata(:)) max(ydata(:))];
    
    if ( isa(h, 'chebtech') )
        % PLOT3(F, G, H)
        
        % Grid data for h:
        data.hGrid.xLine = h.chebpts(npts);
        % Use the maximum of the lenghts of f, g and h to match the number of
        % values returned:
        data.hGrid.xPoints = h.chebpts(len);
        
         % Values on oversampled Chebyshev grid:
        data.zLine = get(prolong(h, npts), 'values');
        data.zPoints = get(prolong(h, len), 'values');
        
    end
    
else

    error('CHEBFUN:CHEBTECH:plotData:dataType', 'Invalid data types.');
    
end

end

function data = plotData(f, g, h)
%PLOTDATA   Useful data values for plotting a SINGFUN object.
%   DATA = PLOTDATA(F) extracts PLOTDATA of the smooth part of F and then scales
%   it by the singular factors given in the EXPONENTS of F.
%
%   DATA = PLOTDATA(F, G) is similar but for plot calls of the form PLOT(F, G),
%   where both F and G are SINGFUN objects.
%
%   DATA = PLOTDATA(F, G, H) is for plots of the form PLOT3(F, G, H). In this
%   instance, DATA also contains fields zLine and zPoints for the data
%   corresponding to H.
%
% See also PLOT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%
% Note: It is assumed that F, G and H are SINGFUNs.

% Make the appropriate call to plotData for the SMOOTHPARTs:
if ( nargin == 1 )
    
    g = [];
    h = [];
    data = plotData(f.smoothPart);
    
elseif (nargin == 2)
    
    h = [];
    f = singfun(f);
    g = singfun(g);
    data = plotData(f.smoothPart, g.smoothPart);
    
elseif (nargin == 3)
    
    f = singfun(f);
    g = singfun(g);
    h = singfun(h);
    data = plotData(f.smoothPart, g.smoothPart, h.smoothPart);
    
end

scaleData = @(x, y, exps) y .* (1+x).^exps(1) .* (1-x).^exps(2);

if ( isempty(g) )
    % PLOT(F):
    
    % Scale the y-data:
    data.yLine = scaleData(data.xLine, data.yLine, f.exponents);
    % Scale the sample point y-data:
    data.yPoints = scaleData(data.xPoints, data.yPoints, f.exponents);
end
    
if ( isa(g, 'singfun') )
    % PLOT(F, G)
    
    % Acquire the grid data of f and scale appropriately:
    data.xLine = scaleData(data.fGrid.xLine, data.xLine, f.exponents);
    data.xPoints = scaleData(data.fGrid.xPoints, data.xPoints, f.exponents);
    
    % Acquire the grid data of g and scale appropriately:
    data.yLine = scaleData(data.gGrid.xLine, data.yLine, g.exponents);
    data.yPoints = scaleData(data.gGrid.xPoints, data.yPoints, g.exponents);
end

if ( isa(h, 'singfun') )
    % PLOT3(F, G, H)

    % Acquire the grid data of h and scale appropriately:
    data.zLine = scaleData(data.hGrid.xLine, data.zLine, h.exponents);
    data.zPoints = scaleData(data.hGrid.xPoints, data.zPoints, h.exponents);
end

% Set the y-limits to something sensible:
data.yLim = getYLimits(data.yLine, f.exponents);

% If F is blowing up, do not use the yLim chosen by Matlab built-in plot:
if ( any(f.exponents < 0) )
    data.defaultYLim = 0;
end

end

function yLim = getYLimits(vals, exps)
%GETYLIMITS   Select y-limits for SINGFUN plots based upon standard deviation.

vals(~isfinite(vals)) = [];
len = length(vals);
mask = true(size(vals));

% If the function blows up on the left, ignore values on the left which are too
% close to the singularity and which will skew the y-limits to be too large:
if ( exps(1) < 0 )
    scl = min(-.2*exps(1), .5);
    numChuck = max(ceil(scl*len), 5);
    mask(1:numChuck) = false;
end

% Same thing on the right:
if ( exps(2) < 0 )
    scl = min(-.2*exps(2), .5);
    numChuck = max(ceil(scl*len), 5);
    mask(end-numChuck:end) = false;
end

% Adjust the y-limits:
masked = vals(mask);
sd = std(masked);
bot = max(min(vals), min(masked) - sd);
top = min(max(vals), max(masked) + sd);
yLim = [bot, top];

end

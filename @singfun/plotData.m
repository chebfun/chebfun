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

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%
% Note: It is assumed that F, G and H are SINGFUNs.

% Make the appropriate call to plotData for the SMOOTHPARTs:
if ( nargin == 1 )
    g = []; h = [];
    data = plotData(f.smoothPart);
elseif (nargin == 2)
    h = [];
    data = plotData(f.smoothPart, g.smoothPart);
elseif (nargin == 3)
    data = plotData(f.smoothPart, g.smoothPart, h.smoothPart);
end

unbnd = false;
scaleData = @(x, y, exps) y .* (1+x).^exps(1) .* (1-x).^exps(2);

if ( isempty(g) )
    % PLOT(F):
    
    % Scale the y-data:
    data.yLine = scaleData(data.xLine, data.yLine, f.exponents);
    % Scale the sample point y-data:
    data.yPoints = scaleData(data.xPoints, data.yPoints, f.exponents);
    
    if ( any(f.exponents < 0 ) )
        unbnd = true;
    end
    
elseif ( isa(g, 'singfun') )
    % PLOT(F, G)
    
    % Acquire the grid data of f and scale appropriately:
    data.xLine = scaleData(data.fGrid.xLine, data.xLine, f.exponents);
    data.xPoints = scaleData(data.fGrid.xPoints, data.xPoints, f.exponents);
    
    % Acquire the grid data of g and scale appropriately:
    data.yLine = scaleData(data.gGrid.xLine, data.yLine, g.exponents);
    data.yPoints = scaleData(data.gGrid.xPoints, data.yPoints, g.exponents);
    
    if ( isa(h, 'singfun') )
        % PLOT3(F, G, H)
        
        % Acquire the grid data of h and scale appropriately:
        data.zLine = scaleData(data.hGrid.xLine, data.zLine, h.exponents);
        data.zPoints = scaleData(data.hGrid.xPoints, data.zPoints, h.exponents);
        
    end
    
    if ( any(g.exponents < 0 ) )
        unbnd = true;
    end
    
end

if ( unbnd )
    % Auto adjust y limits based upon standard deviation
    
    % TODO: This needs much more work!
    
    gl = data.yLine;
    gl(~isfinite(gl)) = [];
    exps = f.exponents;
    l = length(gl);
    mask = true(size(gl));
    scl = min(-.2*exps(1), .5);
    numChuck = max(ceil(scl*l), 5);
    if ( exps(1) >= 0 )
        numChuck = 0;
    end
    mask(1:numChuck) = false;
    scl = min(-.2*exps(2), .5);
    numChuck = max(ceil(scl*l), 5);
    if ( exps(2) >= 0 )
        numChuck = -1;
    end
    mask(end-numChuck:end) = false;
    masked = gl(mask);
    sd = std(masked);
    bot = max(min(gl), min(masked) - sd);
    top = min(max(gl), max(masked) + sd);
    data.yLim = [bot, top];
%     data.yLine(gl > 1.1*top) = NaN;
%     data.yLine(gl < 1.1*bot) = NaN;
end

end
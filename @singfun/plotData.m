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
    g = [];
    data = plotData(f.smoothPart);
elseif (nargin == 2)
    h = [];
    data = plotData(f.smoothPart, g.smoothPart);
elseif (nargin == 3)
    data = plotData(f.smoothPart, g.smoothPart, h.smoothPart);
end

if ( isempty(g) )       
    % PLOT(F):
    % Scale the y-data:
    data.yLine = data.yLine.*(1 + data.xLine).^f.exponents(1).*(1 - data.xLine).^f.exponents(2);
    % Scale the sample point y-data:
    data.yPoints = data.yPoints.*(1 + data.xPoints).^f.exponents(1).*(1 - data.xPoints).^f.exponents(2);
elseif ( isa(g, 'singfun') )   
    % PLOT(F, G)
    
    % Acquire the grid data of f and scale appropriately:
    xLine = data.fGrid.xLine;
    xPoints = data.fGrid.xPoints;
    data.xLine = data.xLine.*(1 + xLine).^f.exponents(1).*(1 - xLine).^f.exponents(2);
    data.xPoints = data.xPoints.*(1 + xPoints).^f.exponents(1).*(1 - xPoints).^f.exponents(2);
           
    % Acquire the grid data of g and scale appropriately:
    xLine = data.gGrid.xLine;
    xPoints = data.gGrid.xPoints;
    data.yLine = data.yLine.*(1 + xLine).^g.exponents(1).*(1 - xLine).^g.exponents(2);
    data.yPoints = data.yPoints.*(1 + xPoints).^g.exponents(1).*(1 - xPoints).^g.exponents(2);    
    
    if ( isa(h, 'singfun') )
        % PLOT3(F, G, H)
        
        % Acquire the grid data of h and scale appropriately:
        xLine = data.hGrid.xLine;
        xPoints = data.hGrid.xPoints;
        data.zLine = data.zLine.*(1 + xLine).^g.exponents(1).*(1 - xLine).^g.exponents(2);
        data.zPoints = data.zPoints.*(1 + xPoints).^g.exponents(1).*(1 - xPoints).^g.exponents(2);
    end
end

end
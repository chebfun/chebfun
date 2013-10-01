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
if ( nargin == 1 )
    g = [];
    data = plotData(f.smoothPart);
elseif (nargin == 2)
    data = plotData(f.smoothPart, g.smoothPart);
elseif (nargin == 3)
    data = plotData(f.smoothPart, g.smoothPart, h.smoothPart);
end

if ( isempty(g) )       
    % PLOT(F):
    % Update the y-data:
    data.yLine = data.yLine.*(data.xLine + 1).^f.exponents(1).*(1 - data.xLine).^f.exponents(2);
    % Update sample point y-data:
    data.yPoints = data.yPoints.*(data.xPoints + 1).^f.exponents(1).*(1 - data.xPoints).^f.exponents(2);
elseif ( isa(g, 'singfun') )   
    % PLOT(F, G)
    
   
    
    if ( isa(h, 'singfun') )
        % PLOT3(F, G, H)
   
    
    
    else
        error('CHEBFUN:SINGFUN:plotdata:DataType', ...
            'Invalid data types.');
        
    end
    
else
    error('CHEBFUN:SINGFUN:plotdata:DataType', ...
        'Invalid data types.');    
end
% Scale data according to exponents.
data = scaleData(data);
% Update the y-data:
x = data.xLine;
y = data.yLine;
z = data.zLine;

y = y.*(x + 1).^f.exponents(1);
y = y.*(1 - x).^f.exponents(2);
data.yLine = y;

% Update sample point y-data:
y = data.yPoints;
x = data.xPoints;
y = y.*(x + 1).^f.exponents(1);
y = y.*(1 - x).^f.exponents(2);
data.yPoints = y;



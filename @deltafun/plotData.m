function data = plotData(f, g)
%PLOTDATA   Useful data values for plotting a DELTAFUN object.
%   DATA = PLOTDATA(F) extracts PLOTDATA of the funPart of F
%   and then appends to it by the data used for delta function plotting.
%
% See also FUN/PLOTDATA.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%
if ( nargin == 1 )
    g = [];
end

% Get plot data from the smooth parth
if ( isempty(g) )
    data = plotData(f.funPart);
else
    data = plotData(f.funPart, g.funPart);
end

% [TODO]: This is commented out since chebfun plot will figure out about
% deltafunctions
%
% % Append vertical lines for deltafunctions
% if ( ~isempty(f.impulses) )
%     % We can only plot the delta functions, not their derivatives:
%     deltaMag = f.impulses(1, :).';
%     deltaLoc = f.location.';
%     
%     % Make vertical line data:
%     deltaY = zeros(2*length(deltaMag), 1);
%     deltaY(1:2:end) = zeros(length(deltaMag), 1);
%     deltaY(2:2:end) = deltaMag;
%     deltaX = zeros(2*length(deltaMag), 1);
%     deltaX(1:2:end) = deltaLoc;
%     deltaX(2:2:end) = deltaLoc;
%     
%     % Append it to existing line data:
%     data.xLine = [data.xLine; deltaX];
%     data.yLine = [data.yLine; deltaY];
%     
%     % Sort the data
%     [data.xLine, index] = sort(data.xLine); 
%     data.yLine = data.yLine(index);
%     
% end
    
end
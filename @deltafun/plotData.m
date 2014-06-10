function data = plotData(f, g)
%PLOTDATA   Useful data values for plotting a DELTAFUN object.
%   DATA = PLOTDATA(F) extracts PLOTDATA of the funPart of F
%   and then appends to it by the data used for delta function plotting.
%
% See also FUN/PLOTDATA.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%
if ( nargin == 1 )
    data = plotData(f.funPart);
    [xDelta, yDelta] = getDeltaData(f);    
else
    % PLOT(F, G)
    % Make sure f has no delta functions.
    if ( isa(f, 'deltafun') )        
        f = f.funPart;
    end    
    
    if ( isa(g, 'deltafun') )
        data = plotData(f, g.funPart);
    else
        data = plotData(f, g);        
    end
    [xDelta, yDelta] = getDeltaData(g);        
end

% Update data with delta functions:
data.xDeltas = xDelta;
data.yDeltas = yDelta;
end

function [xData, yData] = getDeltaData(f)
% GETDELTADATA   Extract data for delta function plotting.

% Initialize empty data:
xData = [];
yData = [];

% Trivial case:
if ( ~isa(f, 'deltafun') )
    return;
end

% Handle delta functions (Derivatives of Delta-functions are not plotted):
if ( ~isempty(f.deltaLoc) )
    % Remove higher derivatives of delta-functions from f:
    f.deltaMag = f.deltaMag(1, :);
    f = simplifyDeltas(f);
    if ( ~isa(f ,'deltafun') ) 
        % No zeroth order delta functions, return:
        return;
    else
        % There are delta functions, prepare data for plotting:                
        xData = f.deltaLoc.';
        % f.deltaMag is necessarily a row vector by now.
        yData = f.deltaMag.';
    end
end
    
end

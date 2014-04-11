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

% Get plot data from the funParth
if ( isempty(g) )
    data = plotData(f.funPart);
else
    data = plotData(f.funPart, g.funPart);
end

% Handle delta functions (Derivatives of Delta-functions are not plotted):
if ( ~isempty(f.deltaLoc) )
    % Remove higher derivatives of delta-functions from f:
    f.deltaMag = f.deltaMag(1, :);
    f = simplifyDeltas(f);
    deltaLoc = f.deltaLoc;
    deltaMag = f.deltaMag;

    data.xDeltas = zeros(3*length(deltaLoc), 1);
    data.xDeltas(1:3:end) = deltaLoc;
    data.xDeltas(2:3:end) = deltaLoc;
    data.xDeltas(3:3:end) = NaN;
    
    data.yDeltas = zeros(3*length(deltaLoc), 1);
    data.yDeltas(2:3:end) = deltaMag(1, :);
    data.yDeltas(3:3:end) = NaN;
else
    data.xDelta = [];
    data.yDelta = [];
end
    
end
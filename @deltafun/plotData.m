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
    g = [];
end

% Get plot data from the funParth
if ( isempty(g) )
    data = plotData(f.funPart);
else
    data = plotData(f.funPart, g.funPart);
end

% Initialise fields for holding data:
data.xDeltas = [];
data.yDeltas = [];

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
        deltaLoc = f.deltaLoc;
        deltaMag = f.deltaMag;
        
        data.xDeltas = zeros(length(deltaLoc), 1);
        data.xDeltas = deltaLoc.';
        
        data.yDeltas = zeros(length(deltaLoc), 1);
        data.yDeltas = deltaMag(1, :).';
    end
end
    
end

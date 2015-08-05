function data = dispData(f)
%DISPDATA   Useful DELTAFUN information for DISPLAY at higher levels.
%   DATA = DISPDATA(F) extracts useful information from the given DELTAFUN F and
%   the information DATA will be used by DISPLAY at higher levels.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

data = dispData(f.funPart);

if ( anyDelta(f) )    
    deltaData.name = 'deltas';
    deltaData.data.deltaMag = f.deltaMag;
    deltaData.data.deltaLoc = f.deltaLoc;
    data = [data {deltaData}];
end

end

function data = dispData(f)
%DISPDATA   Useful DELTAFUN information for DISPLAY at higher levels.
%   DATA = DISPDATA(F) extracts useful information from the given DELTAFUN F and
%   the information DATA will be used by DISPLAY at higher levels. Currently,
%   the only information it extracts is exponents.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( anyDelta(f) )    
    data{1}.name = 'deltas';
    data{1}.data.deltaMag = f.deltaMag;
    data{1}.data.deltaLoc = f.deltaLoc;
else
    data{1} = [];
end

end
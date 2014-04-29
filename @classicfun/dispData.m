function data = dispData(f)
%DISPDATA   Useful information for DISPLAY at higher levels.
%   DATA = DISPDATA(F) extracts useful information from the given CLASSICFUN 
%   F and the information will be used by DISPLAY at higher levels. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

data = dispData(f.onefun);

% More information for F can be appended to DATA:

end

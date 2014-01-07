function info = dispInfo(f)
%DISPINFO   Useful information for DISPLAY at higher levels.
%   INFO = DISPINFO(F) extracts useful information from the given FUN F and
%   the information INFO will be used by DISPLAY at higher levels. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

info = dispInfo(f.onefun);

% More information for F can be appended to INFO:

end
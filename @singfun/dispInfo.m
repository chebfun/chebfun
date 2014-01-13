function info = dispInfo(f)
%DISPINFO   Useful information for DISPLAY at higher levels.
%   INFO = DISPINFO(F) extracts useful information from the given SINGFUN F and
%   the information INFO will be used by DISPLAY at higher levels. Currently, 
%   the only information it extracts is exponents.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

info{1}.name = 'exponents';
info{1}.data = f.exponents;

% More information can be appended to INFO:

end
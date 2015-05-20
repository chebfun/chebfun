function [minVal, minPos] = min(f)
%MIN   Global minimum of a DELTAFUN.
%   MINVAL = MIN(F) returns the global minimum of the DELTAFUN F.
%
%   [MINVAL, MINPOS] = MIN(F) returns also a value such that MINVAL = F(MINPOS).
%
% See also MAX, MINANDMAX.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% To avoid code duplication, we simply call MINANDMAX():
[minVal, minPos] = minandmax(f);
minVal = minVal(1,:);
minPos = minPos(1,:);

end

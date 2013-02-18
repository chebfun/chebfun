function [minVal, minPos] = min(f)
%MIN	Global minimum of a FUNCHEB2 on [-1,1].
%   MINVAL = MIN(F) is the global minimum MINVAL of the FUNCHEB2 F on [-1,1].
%
%   [MINVAL, MINPOS] = MIN(F) returns also the value MAXPOS such that MINVAL =
%   G(MINPOS).
%
% See also MAX, MINANDMAX.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% To avoid code duplication, we simply call minandmax:
[minVal, minPos] = minandmax(f);
minVal = minVal(1,:);
minPos = minPos(1,:);

end
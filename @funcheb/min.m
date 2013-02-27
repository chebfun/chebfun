function [minVal, minPos] = min(f)
%MIN	Global minimum of a FUNCHEB on [-1,1].
%   MINVAL = MIN(F) is the global minimum MINVAL of the FUNCHEB F on [-1,1].
%
%   [MINVAL, MINPOS] = MIN(F) returns also the value MAXPOS such that MINVAL =
%   F(MINPOS).
%
%   If F is complex-valued then absolute values are taken to determine maxima
%   but the resulting value corresponds to that of the original function. That
%   is MINVAL = feval(F, MINPOS) where [~, MINPOS] = MIN(abs(F)).
%
% See also MAX, MINANDMAX.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% To avoid code duplication, we simply call minandmax:
[minVal, minPos] = minandmax(f);
minVal = minVal(1,:);
minPos = minPos(1,:);

end
function [maxVal, maxPos] = max(f)
%MAX	Global maximum of a FUNCHEB2 on [-1,1].
%   MAXVAL = MAX(F) is the global maximum MAXVAL of the FUNCHEB2 F on [-1,1].
%
%   [MAXVAL, MAXPOS] = MAX(F) returns also the value MAXPOS such that MAXVAL =
%   G(MAXPOS).
%
% See also MIN, MINANDMAX.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% To avoid code duplication, we simply call minandmax:
[maxVal, maxPos] = minandmax(f);
maxVal = maxVal(2,:);
maxPos = maxPos(2,:);

end


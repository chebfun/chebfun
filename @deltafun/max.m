function [maxVal, maxPos] = max(f)
%MAX   Global maximum of a DELTAFUN.
%   MAXVAL = MAX(F) returns the global maximum of the DELTAFUN F.
%
%   [MAXVAL, MAXPOS] = MAX(F) returns also a value such that MAXVAL = F(MAXPOS).
%
% See also MIN, MINANDMAX.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% To avoid code duplication, we simply call MINANDMAX():
[maxVal, maxPos] = minandmax(f);
maxVal = maxVal(2,:);
maxPos = maxPos(2,:);

end

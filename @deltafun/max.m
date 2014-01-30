function [maxVal, maxPos] = max(f)
%MAX   Global maximum of a DELTAFUN on [-1,1].
%   MAXVAL = MAX(F) returns the global maximum of the DELTAFUN F on [-1,1].
%
%   [MAXVAL, MAXPOS] = MAX(F) returns also a value such that MAXVAL = F(MAXPOS).
%
%   [TODO]: If F is complex-valued then absolute values are taken to determine
%   maxima but the resulting value corresponds to that of the original function.
%   That is, MAXVAL = feval(F, MAXPOS) where [~, MAXPOS] = MAX(abs(F));
%
% See also MIN, MINANDMAX.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% To avoid code duplication, we simply call MINANDMAX():
[maxVal, maxPos] = minandmax(f);
maxVal = maxVal(2,:);
maxPos = maxPos(2,:);

% TODO: NH: DOCS are for SINGFUN!

end

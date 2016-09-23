function [maxVal, maxPos] = max(f)
%MAX   Global maximum of a CLASSICFUN on [a,b].
%   MAXVAL = MAX(F) returns the global maximum of the CLASSICFUN F on [a,b].  If F is
%   an array-valued CLASSICFUN, MAXVAL is a row vector whose Kth entry is the global
%   maximum of the Kth column of F.
%
%   [MAXVAL, MAXPOS] = MAX(F) returns also a value such that MAXVAL = F(MAXPOS).
%
%   If F is complex-valued then absolute values are taken to determine maxima
%   but the output value corresponds to that of the original function. That
%   is, MAXVAL = feval(F, MAXPOS) where [~, MAXPOS] = MAX(abs(F));
%
% See also MIN, MINANDMAX.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% To avoid code duplication, we simply call MINANDMAX:
[maxVal, maxPos] = minandmax(f);
maxVal = maxVal(2,:);
maxPos = maxPos(2,:);

end

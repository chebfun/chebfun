function f = thresholdBreakpointValues(f)
%THRESHOLDBREAKPOINTVALUES   Set small breakpoint values to zero.
%   G = THRESHOLDBREAKPOINTVALUES(F), where F is a CHEBFUN, returns a CHEBFUN
%   G such that all breakpoint values smaller than VSCALE(F)*EPSLEVEL(F) are
%   set to zero.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

breakVals = f.pointValues;
breakVals(abs(breakVals) < vscale(f)*epslevel(f)) = 0;
f.pointValues = breakVals;

end

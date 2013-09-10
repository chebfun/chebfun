function f = introduceBreakpoint(f, x)
%INTRODUCEBREAKPOINT    Introduce a breakpoint in the domain of a CHEBFUN.
%   INTRODUCEBREAKPOINT(F, X) introduces breakpoints in the CHEBFUN F at the
%   points in X.
%
% See also DEFINEPOINT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% This function is simply a wrapper for @CHEBFUN/DEFINEPOINT().
x = x(:);
f = definePoint(f, x, feval(f, x));

end

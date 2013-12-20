function val = feval(f, x)
%FEVAL   Evaluate a DELTAFUN.
%   FEVAL(F, X) evaluates the DELTAFUN F at the given points X.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% For evaluation, the underlying smoothfun, i.e. the smoothPart is evaluated
% at X first and then the values computed are scaled by the singular factors.

% Evaluate the smooth part.
val = feval(f.funPart, x);

% [TODO]: finish this:
end

function val = feval(f, x)
%FEVAL   Evaluate a SINGFUN.
%   FEVAL(F, X) evaluates the SINGFUN F at the given points X.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% For evaluation, the underlying SMOOTHFUN, i.e., f.smoothPart, is first
% evaluated at X and then the values computed are scaled by the singular
% factors.

% Evaluate the smooth part:
val = feval(f.smoothPart, x);

% Multiply now with the singular factors:
if ( f.exponents(1) )
    % If there is a non-trivial left singular factor:
    val = val.*(1 + x).^f.exponents(1);
end

if ( f.exponents(2) )
    % If there is a non-trivial right singular factor:
    val = val.*(1 - x).^f.exponents(2);
end

end

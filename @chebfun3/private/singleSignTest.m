function [ss, wzero, ispos] = singleSignTest(F)
%SINGLESIGNTEST   Heuristic check for sign change.
%   SINGLESIGNTEST(F) returns 1 if the values of F on a tensor grid are of 
%   the same sign.
%
%   [SS, WZERO, ISPOS] = SINGLESIGNTEST(F) returns WZERO = 1 if a zero has 
%   been found and ISPOS = 1 if F is positive-valued.
%
%   The algorithm works by sampling F on a tensor grid and checking if
%   those values are of the same sign. This command is mainly for internal 
%   use in CHEBFUN3 commands.
%
% See also CHEBFUN3/ABS, CHEBFUN3/SQRT and CHEBFUN3/LOG.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

tol = 2*chebfun3eps;
ss = false;                  % Assume false
ispos = false;

X = sample(F);  % Evaluate on a grid using FFTs.
X = X(:);

if ( all( X > -tol * F.vscale ))  % If all values are nonnegative     
    ss = true;
    ispos = true;
elseif ( all( X < tol * F.vscale))        % If all values are not positive
    ss = true;
end

wzero = any(X == 0);        % Any exact zeros on the grid?

end
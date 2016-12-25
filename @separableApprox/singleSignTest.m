function [out, wzero] = singleSignTest(F)
%SINGLESIGNTEST   Heuristic check to see if F changes sign.
% 
%   SINGLESIGNTEST(F) returns 1 if the values of F on a tensor grid are of the
%   same sign.
%
%   [OUT, WZERO] = SINGLESIGNTEST(F), returns WZERO = 1 if a zero has been
%   found.
%
%   The algorithm works by sampling F on a tensor-grid and checking if
%   those values are of the same sign. This command is mainly for internal use
%   in SEPARABLEAPPROX commands.
%
% See also ABS, SQRT, LOG. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

tol = chebfun2eps;
out = false;                  % Assume false

% Evaluate on a grid:
X = sample(F);

X = X(:);

if ( all( X >= -tol * F.vscale ))   % If all values are nonnegative         
    out = true;  
elseif ( all( X <= tol * F.vscale))  % If all values are not positive
    out = true; 
end

wzero = any(X == 0);        % Any exact zeros on the grid?

end
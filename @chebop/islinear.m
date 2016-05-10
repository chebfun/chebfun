function out = islinear(N, u)
%ISLINEAR   Determine linearity of a CHEBOP.
%   OUT = ISLINEAR(N) determines the linearity of the CHEBOP N. In particular:
%       OUT(1) = 1 if N.OP is linear, 0 otherwise.
%       OUT(2) = 1 if N.LBC is linear, 0 otherwise.
%       OUT(3) = 1 if N.RBC is linear, 0 otherwise.
%       OUT(4) = 1 if N.BC is linear, 0 otherwise.
%
%   OUT = ISLINEAR(N, U) performs the linearization of N around the function U,
%   rather than the zero function.
%
% See also LINEARIZE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

if ( nargin < 2 )
    % Allow LINEARIZE() to chose the linearization variable:
    u = [];
end

% Call LINEARIZE(). It might occur that something non-differentiable appears in
% the .op field of the chebop, e.g. @(u) diff(u,2) + abs(u). This will cause
% an error down the stack, since abs is not Frechet differentiable, and the
% corresponding ADCHEBFUN method throws an error. Thus, we wrap the call to
% linearize in a try-catch statement, and listen out for the 'nonDifferentiable'
% error identifier.
try 
    [ignored1, ignored2, out] = linearize(N, u, [], 0);
catch ME
    if ( ~isempty(strfind(ME.identifier, 'notDifferentiable')) || ...
        strcmp(ME.identifier, 'CHEBFUN:CHEBOP:linearize:invalidInitialGuess') )
        % Something non-differentiable appeared in the operator, so it can't be
        % linear! Alternatively, the operator fails to evaluate on the initial
        % guess passed, so it can't in that case be linear either!
        out = false;
    else
        rethrow(ME);
    end
end

end

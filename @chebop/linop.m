function [L, f, FAIL] = linop(N)
%LINOP   Convert a CHEBOP to a LINOP.
%   L = LINOP(N) converts a CHEBOP N to a linop L if N is a linear operator. If
%   N is not linear, an error message is returned.
%
%   [L, F] = LINOP(N) returns also the affine part F of the linear CHEBOP N such
%   that L*u + F(x) = N.op(x,u).
%
%   [L, F, FAIL] = LINOP(N) will prevent an error from being thrown if N is not
%   linear, but instead return FAIL = TRUE and L = []. If N is linear, FAIL =
%   FALSE.
%
% See also LINOP, CHEBOP/LINEARIZE, CHEBOP/ISLINEAR.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Tell CHEBOP/LINEARIZE() to stop if it detects nonlinearity:
linCheck = 1; 

% Linearize, thereby obtaining linearity information and a LINOP:
[L, f, isLinear] = linearize(N, N.init, [], linCheck);

% We need the entire operator (including BCs) to be linear:
FAIL = ~all(isLinear);

% Throw an error is the CHEBOP is nonlinear:
if ( FAIL && (nargout < 3) )
    error('CHEBFUN:CHEBOP:linop:nonlinear',...
        'This does not appear to be a linear operator.')
end

% All values of the LINOPCONSTRAINT stored in L will be of incorrect sign when
% returned from LINEARIZE(), if we want to use it for a LINOP backslash. This is
% because when problems are solved with LINOP backslash, the solution to the
% problem is the output itself, while in a Newton iteration, we have to add the
% output of the LINOP solution to the current guess. Thus, flip the signs of the
% values of L.constraint.
L.constraint = -L.constraint;

end

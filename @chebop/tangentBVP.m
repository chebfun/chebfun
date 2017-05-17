function [t, tau] = tangentBVP(N, u, lambda, told, tauold, x, dSum, prefs)
%TANGENTBVP    Find a tangent to the solution curve at a given point
%
% TANGENTBVP finds a tangent to the curve 
%   H(U, LAMBDA) = 0
% at the given point, by solving a linear boundary-value problem.
%
% Calling sequence:
%   [T, TAU] = TANGENTBVP(N, U, LAMBDA, TOLD, TAUOLD, X, DSUM, PREFS)
%
% Here, the inputs are:
%   
%   N       : A chebop, whose N.op arguments are x, u and lambda, and boundary
%             conditions also depend on u and lambda.
%   U       : The current solution on the solution curve.
%   LAMBDA  : The value of lambda for the solution U on the curve.
%   TOLD    : Previous tangent function.
%   TAUOLD  : Previous tangent scalar.
%   X       : The independent CHEBFUN on the domain of the problem.
%   DSUM    : A diagonalised sum operator on the domain of the problem.
%   PREFS   : A CHEBOPPREF object.
%
% The output is a tangent to the solution curve, which is an "Inf + 1"
% dimensional object. It consists of:
%   T   : A CHEBFUN.
%   TAU : A scalar.
%
% See also: followpath, newtonBVP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% The domain of the problem
d = domain(u);

% Linearize the augmented operator around the current solution:
S = linearize(N, [u; lambda], x);

% Create a chebmatrix for the constraint which goes at the bottom of the linear
% operator:
Jm = [told'*dSum, tauold];

% Store the constraint of S (as it will be cleared out below when we recreate
% the linop):
Scon = S.constraint;

% Derivative of augmented operator with the constraint at the bottom:
S = linop([S;Jm]);

% Reassign constraint. For the tangent problem, we always have homogenous
% conditions:
Scon.values = 0*Scon.values;
S.constraint = Scon;

% Right hand side
rhs = [chebfun(0, d); 1];

% Solve for the tangent:
ttau = linsolve(S, rhs, prefs);

% Extract the function and scalar:
t = ttau{1};
tau = ttau{2};

% Normalize the tangent returned:
scale = sqrt(t'*t + tau^2);
t = t/scale;
tau = tau/scale;

end
function [u, info] = solvebvpLinear(L, rhs, pref, displayInfo)
%SOLVEBVP  Solve a linear CHEBOP BVP system.
%
% [U, INFO] = SOLVEBVPLINEAR(L, RHS, PREF, DISPLAYINFO), where:
%   L is a linear CHEBOP
%   RHS is a CHEBMATRIX
%   PREF is a CHEBOPPREF
%   DISPLAYINFO is a function handle
%
% attempts to solve the linear BVP
%
%       L*U = RHS + boundary conditions specified by L
%
% The output U is a CHEBMATRIX, and INFO is a MATLAB struct with useful
% information. This method should generally not called directly by the user,
% which should rather call the CHEBOP/MLDIVIDE or CHEBOP/SOLVEBVP methods.
%
% Observe that when this method is called, any affine parts of the original
% CHEBOP have been absorbed into RHS.
%
% See also: chebop/solvebvp

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get defaults:
if ( nargin < 3 )
    pref = cheboppref();
end

% All values of the LINOPCONSTRAINT stored in L will be of incorrect sign when
% returned from LINEARIZE(), if we want to use it for a LINOP backslash. This is
% because when problems are solved with LINOP backslash, the solution to the
% problem is the output itself, while in a Newton iteration, we have to add the
% output of the LINOP solution to the current guess. Thus, flip the signs of the
% values of L.constraint:
L.constraint = flipSigns(L.constraint);

% Solutions to the linearized problems need to be more accurate than the
% nonlinear iteration tolerance.
linpref = pref;
linpref.errTol = max( eps, pref.errTol/100 );

% Solve the linear problem:
u = linsolve(L, rhs, linpref);

% Norm of residual:
normRes = norm(L*u - rhs, 'fro');

if ( nargin > 4 )
    % Print information after linear problem has been solved:
    displayInfo('linear', u, normRes, pref)
end

if ( nargout > 1 )
    % Return the norm of the residual in the INFO struct:
    info.error = normRes;
end

end
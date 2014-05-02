function [u, info] = solvebvpLinear(L, rhs, res, pref, displayInfo)

% TODO: Document.

% TODO: Pass RHS - RES as RHS?

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

% Solve the linear problem:
u = linsolve(L, rhs - res, pref);

% Norm of residual. Any affine parts will be stored in the RES variable, so need
% to subtract that from RHS to get the correct answer:
resNorm = norm(L*u - (rhs - res), 'fro');

if ( nargin > 4 )
    % Print information after linear problem has been solved:
    displayInfo('linear', u, resNorm, pref)
end

if ( nargout > 1 )
    % Return the norm of the residual in the INFO struct:
    info.error = resNorm;
end

end
function [u, info] = solvebvpLinear(L, rhs, pref, displayInfo)

% TODO: Document.

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
u = linsolve(L, rhs, pref);

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
function [u, info] = solvebvpLinear(N, L, rhs, residual, x, displayInfo, pref)

% All values of the LINOPCONSTRAINT stored in L will be of incorrect sign
% when returned from LINEARIZE(), if we want to use it for a LINOP
% backslash. This is because when problems are solved with LINOP backslash,
% the solution to the problem is the output itself, while in a Newton
% iteration, we have to add the output of the LINOP solution to the current
% guess. Thus, flip the signs of the values of L.constraint:
L.constraint = flipSigns(L.constraint);

% Solve the linear problem
u = linsolve(L, rhs - residual, pref);

% TODO: Return residual as well?
uBlocks = u.blocks;

% Norm of residual. Any affine parts will be stored in the RESIDUAL variable, so
% need to subtract that from RHS to get the correct answer.
normRes = norm(chebfun(L*u - (rhs-residual)), 2);

% Print information after linear problem has been solved
displayInfo('linear', u, normRes, pref)

% Return the norm of the residual in the INFO struct.
info.error = normRes;

end
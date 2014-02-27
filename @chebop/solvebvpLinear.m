function [u, info] = solvebvpLinear(N, L, rhs, residual, x, displayInfo, pref)

% All values of the LINOPCONSTRAINT stored in L will be of incorrect sign
% when returned from LINEARIZE(), if we want to use it for a LINOP
% backslash. This is because when problems are solved with LINOP backslash,
% the solution to the problem is the output itself, while in a Newton
% iteration, we have to add the output of the LINOP solution to the current
% guess. Thus, flip the signs of the values of L.constraint:
L.constraint = flipSigns(L.constraint);
% TODO: Pass in preferences?

% Solve the linear problem
u = linsolve(L, rhs - residual);

% TODO: Return residual as well?
uBlocks = u.blocks;

% Norm of residual
normRes = norm(chebfun(L*u - rhs));

% Print information after linear problem has been solved
displayInfo('linear', u, normRes, pref)

% TODO: Probably want a norm method for chebmatrices. THIS WILL BREAK IN
% CASE OF SYSTEMS.
err = feval(N.op, x , uBlocks{:}) - rhs;
if isa(err, 'chebmatrix')
    err = err.blocks;
    info.error = norm(norm(err{:}));
else
    info.error = norm(err);
end

end
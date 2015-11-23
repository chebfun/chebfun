function [u, lam, iter, retract] = newtonBVP(N, u, lam, t, tau, x, dsum, prefs)
%NEWTONBVP    Get back to the solution curve after taking a tangent step
%
% TANGENTBVP finds a correction after taking a tangent step to get back on the
% curve
%   N(U, LAMBDA) = 0
% from the given point (UDOT, LAMDOT), by solving a nonlinear boundary-value
% problem.
%
% Calling sequence:
%   [U, LAM, ITER, RETRACT] = NEWTONBVP(N, U, LAM, T, TAU, X, DSUM, PREFS)
%
% Here, the inputs are:
%   
%   N       : A chebop, whose N.op arguments are x, u and lambda, and boundary
%             conditions also depend on u and lambda.
%   U       : The current approximate solution (where we ended up after taking
%             the tangent step).
%   LAM     : The current approximate path following parameter (where we ended
%             up after taking the tangent step).
%   T       : Previous tangent function (which we want to be orthogonal to).
%   TAU     : Previous tangent scalar (which we want to be orthogonal to).
%   X       : The independent CHEBFUN on the domain of the problem.
%   DSUM    : A diagonalised sum operator on the domain of the problem.
%   PREFS   : A CHEBOPPREF object.
%
% The output, if the algorithm converges, is a point (U, LAM) that lies on the
% solution curve N(U, LAM) = 0. The output order is as follows:
%   U       : A CHEBFUN that lies on the solution curve N(U, LAM) = 0 for the
%             LAM value returned.
%   LAM     : A scalar that lies on the solution curve N(U, LAM) = 0 for the
%             U value returned.
%   ITER    : The number of Newton iteration required to converge.
%   RETRACT : A Boolean, that tells the followpath algorithm whether it should
%             retract (take a smaller tangent step).
%
% See also: followpath, newtonBVP.

% Begin by creating a chebmatrix for the constraint which goes at the
% bottom of the linearized operator during the Newton iteration.
Jm = [t'*dsum, tau];

% Store the norm of current solution to be able to check whether we've
% converged.
normu = norm(u);

% Compute Newton correction
% Do we accept the tangent step?
accept = 0;
% Count number of Newton iterations:
iter = 0;
% Does the Newton correction algorithm want the pathfollowing one to retract
% (take a smaller tangent step)?
retract = 0;

% Do a Newton iteration to get us back on the solution curve.
while ( ~accept )
    
    % Linearise the operator N, obtaining the residual at the same time:
    [S, res] = linearize(N, [u; lam], x);
    
    % Store the current constraints (boundary conditions) as they'll be cleared
    % out when we augment the operator below:
    Scon = S.constraint;
    Scon.values = 0*Scon.values;
    
    % Create a linop by adding the Jm constraint to S:
    L = linop([S; Jm]);
    
    % Reassign the (boundary) constraints:
    L.constraint = Scon;
        
    % RHS is the residual of the differential equation, we then add a 0 at the
    % bottom for the functional condition:
    rhs = [-res{1}; 0];
    
    % Solve the linear system to obtain a Newton correction:
    dudlam = linsolve(L, rhs, prefs);
    
    % Extract the function and scalar parts:
    du = dudlam{1};
    dlam = dudlam{2};
    
    % Take a Newton step. Since we should be close to a solution on the curve,
    % we take a full Newton step.
    u = u + du;
    lam = lam + dlam;

    % Increase the iteration counter.
    iter = iter + 1;
    
    if ( norm(du, 2)/normu < 1e-3 )
        % Hoorayh. We've converged, so accept the tangent step!
        accept = 1;
    elseif ( iter >= 5 )
        % We wanted to take too many iterations, so tell the followpath
        % algorithm to retract and take a smaller tangent step.
        retract = 1;
        return
    end
end
end
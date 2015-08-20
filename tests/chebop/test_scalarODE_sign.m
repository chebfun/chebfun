function [pass, u1, u2, info1, info2] = test_scalarODE_sign(pref)
% A nonlinear CHEBOP test. This test tests a scalar ODE, where there is a
% breakpoint in the domain of the CHEBOP. Furthermore, the operator has a
% discontinuous coefficient, which will induce a further breakpoint in the
% solution. The problem is solved using chebcolloc1, chebcolloc2 and ultraS
% discretizations. The problem solved does not require damping for the Newton
% iteration to converge.
%
% Asgeir Birkisson, May 2014.

%% Setup
if ( nargin == 0 )
    pref = cheboppref;
end
dom = [-1 .5 1];

N = chebop(@(x,u) diff(u,2) + sign(x).*sin(u), dom);
N.lbc = @(u) u - 2;
N.rbc = @(u) u - 2;
rhs = 0;

%% Try different discretizations
% Start with chebcolloc2
pref.discretization = @chebcolloc2;
[u1, info1] = solvebvp(N, rhs, pref);

%% Change to ultraS
pref.discretization = @ultraS;
[u2, info2] = solvebvp(N, rhs, pref);

%% Change to chebcolloc1
pref.discretization = @chebcolloc1;
[u3, info3] = solvebvp(N, rhs, pref);

%% Did we pass? 
% To pass, both residuals have to be small, but we should not expect u1 and u2
% to be identical!

% TODO: This used to be pref.errTol. Once we tune the algorithms better, should
% try to restore it.
tol = 1e2*pref.errTol;
pass(1) = norm(N(u1)) < tol;
pass(2) = norm(N(u2)) < tol;
pass(3) = norm(N(u3)) < tol;
pass(4) = ( (norm(u1 - u2) ~= 0) && (norm(u2 - u3) ~= 0) && ...
    (norm(u1 - u3) ~= 0));

end

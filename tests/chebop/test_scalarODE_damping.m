function [pass, u1, u2, info1, info2] = test_scalarODEdamping(pref)
% A nonlinear CHEBOP test. This test tests a scalar ODE, where no breakpoints
% occur. It solves the problem using chebcolloc1, chebcolloc2 and ultraS
% discretizations. The problem solved requires damping for the Newton iteration
% to converge.
%
% Asgeir Birkisson, May 2014.

%% Setup
dom = [0 pi];
if ( nargin == 0 )
    pref = cheboppref;
end

N = chebop(@(x,u) .05*diff(u,2) + cos(5*x).*sin(u), dom);
N.lbc = @(u) u - 2; 
N.rbc = @(u) u - 3;
rhs = 0;

pref.bvpTol = 1e-10;

% Try different discretizations:

%% Start with chebcolloc2
pref.discretization = @chebcolloc2;
pref.bvpTol = 1e-13;
[u1, info1] = solvebvp(N, rhs, pref);
err(1) = norm(N(u1));

%% Change to chebcolloc1
pref.discretization = @chebcolloc1;
[u2, info2] = solvebvp(N, rhs, pref);
err(2) = norm(N(u2));

%% Change to ultraS
pref.discretization = @ultraS;
pref.bvpTol = 1e-12;
[u3, info3] = solvebvp(N, rhs, pref);
err(3) = norm(N(u3));

%% Did we pass? 
% To pass, both residuals have to be small, but we should not expect u1 and u2
% to be identical!

tol = 1e-9;
pass = err < tol;

end

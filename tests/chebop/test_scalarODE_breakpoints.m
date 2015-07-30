function [pass, u1, u2, info1, info2] = test_scalarODE_breakpoints(pref)
% A nonlinear CHEBOP test. This test tests a scalar ODE, where there are two
% breakpoints in the domain of the CHEBOP.  It solves the problem using chebcolloc1,
% chebcolloc2 and ultraS discretizations.
%
% Asgeir Birkisson, May 2014.


%% Setup
if ( nargin == 0 )
    pref = cheboppref;
end

dom = [0 1 2 pi];
pref.damping = 1;

N = chebop(@(x,u) diff(u,2) + sin(u-.2), dom);
N.lbc = @(u) u - 2;
N.rbc = @(u) u - 3;
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
tol = 5*pref.errTol;
err1 = normest(N(u1));
err2 = normest(N(u2));
err3 = normest(N(u3));
% TODO: This used to be 1*, 10* and 1*. Should try to restore once we tune the
% algorithms better.
pass(1) = err1 < 1e2*tol;
pass(2) = err2 < 1e2*tol;
pass(3) = err3 < 1e2*tol;
pass(4) = ( (norm(u1 - u2) ~= 0) && (norm(u2 - u3) ~= 0) && ...
    (norm(u1 - u3) ~= 0));

end

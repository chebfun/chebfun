function [pass, u1, u2, info1, info2] = test_scalarODE_breakpoints(pref)

if ( nargin == 0 )
    pref = cheboppref;
end

% Setup
dom = [0 1 2 pi];
pref.plotting = 'on';
pref.display = 'iter';
pref.damped = 1;

N = chebop(@(x,u) diff(u,2) + sin(u-.2), dom);
N.lbc = @(u) u-2;
N.rbc = @(u) u - 3;
rhs = 0;

%% Try different discretizations
% Start with collocation -- no further action required
[u1, info1] = solvebvp(N, rhs, pref);

%% Change to ultraS
pref.discretization = @ultraS;
[u2, info2] = solvebvp(N, rhs, pref);

%% Did we pass? 
% To pass, both residuals have to be small, but we should not expect u1 and u2
% to be identical!
tol = pref.errTol;
pass(1) = normest(N(u1)) < tol;
pass(2) = normest(N(u2)) < tol;
pass(3) = ( norm(u1-u2) ~= 0 );

function pass = test_linearScalarODEs(pref)
% A linear CHEBOP test. This test tests a scalar ODE, both with and without
% breakpoints, as well as with discontinuous coefficients. It solves the
% problems using chebcolloc1, chebcolloc2 and ultraS discretizations.

%% Setup
dom = [0 pi];
if ( nargin == 0 )
    pref = cheboppref;
end
%pref.bvpTol = 1e-11;

%% Simple scalar problem
N = chebop(@(x,u) diff(u,2) + x.*u, dom);
N.lbc = 2; 
N.rbc = 3;

x = chebfun(@(x) x, dom);
rhs = sin(x);

%% Try different discretizations
% Start with chebcolloc2
pref.discretization = @chebcolloc2;
u1 = solvebvp(N, rhs, pref);

%% Change to ultraS
pref.discretization = @ultraS;
u2 = solvebvp(N, rhs, pref);

%% Change to chebcolloc1
pref.discretization = @chebcolloc1;
u3 = solvebvp(N, rhs, pref);

%% Did we pass? 
% To pass, both residuals have to be small, but we should not expect u1 and u2
% to be identical!
tol = 1e3*pref.bvpTol;
pass(1) = norm(N(u1)-rhs) < tol && ( u1(0) - 2 < tol) && ( u1(pi) - 3 < tol);
pass(2) = norm(N(u2)-rhs) < tol && ( u2(0) - 2 < tol) && ( u2(pi) - 3 < tol);
pass(3) = norm(N(u3)-rhs) < tol && ( u3(0) - 2 < tol) && ( u3(pi) - 3 < tol);
pass(4) = ( (norm(u1 - u2) ~= 0) && (norm(u2 - u3) ~= 0) && ...
    (norm(u1 - u3) ~= 0));


%% Problem with breakpoints
dom = [-1 0 pi];
N = chebop(@(x,u) diff(u,2) + cos(x).*u, dom);
N.lbc = 2; 
N.rbc = -1;

x = chebfun(@(x) x, dom);
rhs = sin(x);

%% Try different discretizations
% Start with chebcolloc2
pref.discretization = @chebcolloc2;
u4 = solvebvp(N, rhs, pref);

%% Change to ultraS
pref.discretization = @ultraS;
u5 = solvebvp(N, rhs, pref);

%% Change to chebcolloc1
pref.discretization = @chebcolloc1;
u6 = solvebvp(N, rhs, pref);

%% Did we pass? 
% To pass, both residuals have to be small, but we should not expect u3 and u4
% to be identical!
tol = 1e3*pref.bvpTol;
pass(5) = norm(N(u4)-rhs) < tol && ( u4(-1) - 2 < tol) && ( u4(pi) + 1 < tol);
pass(6) = norm(N(u5)-rhs) < tol && ( u5(-1) - 2 < tol) && ( u5(pi) + 1 < tol);
pass(7) = norm(N(u6)-rhs) < tol && ( u6(-1) - 2 < tol) && ( u6(pi) + 1 < tol);
pass(8) = norm(jump(u4, 0)) < tol && norm(jump(u5, 0)) < tol && ...
    norm(jump(u6, 0)) < tol;
pass(9) = ( (norm(u4 - u5) ~= 0) && (norm(u4 - u6) ~= 0) && ...
    (norm(u5 - u6) ~= 0));

%% Problem with discontinuous coefficients
dom = [-1 0 pi];
N = chebop(@(x,u) diff(u,2) + abs(x-1).*u, dom);
N.lbc = 2; 
N.rbc = -1;

x = chebfun(@(x) x, dom);
rhs = sin(x);

%% Try different discretizations
% Start with collocation
pref.discretization = @chebcolloc2;
u7 = solvebvp(N, rhs, pref);

%% Change to ultraS
pref.discretization = @ultraS;
u8 = solvebvp(N, rhs, pref);

%% Change to chebcolloc1
pref.discretization = @chebcolloc1;
u9 = solvebvp(N, rhs, pref);


%% Did we pass? 
% To pass, both residuals have to be small, but we should not expect u3 and u4
% to be identical!
tol = 1e1*pref.bvpTol;
pass(10) = norm(N(u7)-rhs) < 50*tol && ...
    ( u7(-1) - 2 < tol) && ( u7(pi) + 1 < tol);
pass(11) = norm(N(u8)-rhs) < 10*tol && ...
    ( u8(-1) - 2 < tol) && ( u8(pi) + 1 < tol);
pass(12) = norm(N(u9)-rhs) < 10*tol && ...
    ( u9(-1) - 2 < tol) && ( u9(pi) + 1 < tol);
pass(13) = norm(jump(u7, 0)) < tol && norm(jump(u8, 0)) < tol && ...
    norm(jump(u9, 0)) < tol;
pass(14) = ( norm(u5-u6) ~= 0 );

end

function pass = test_linearScalarODEs(pref)

% Setup
dom = [0 pi];
if ( nargin == 0 )
    pref = cheboppref;
end

%% Simple scalar problem
N = chebop(@(x,u) diff(u,2) + x.*u, dom);
N.lbc = 2; 
N.rbc = 3;

x = chebfun(@(x) x, dom);
rhs = sin(x);

%% Try different discretizations
% Start with collocation -- no further action required
u1 = solvebvp(N, rhs, pref);

%% Change to ultraS
pref.discretization = @ultraS;
u2 = solvebvp(N, rhs, pref);

%% Did we pass? 
% To pass, both residuals have to be small, but we should not expect u1 and u2
% to be identical!
tol = pref.errTol;
pass(1) = norm(N(u1)-rhs) < tol && ( u1(0) - 2 < tol) && ( u1(pi) - 3 < tol);
pass(2) = norm(N(u2)-rhs) < tol && ( u2(0) - 2 < tol) && ( u2(pi) - 3 < tol);
pass(3) = ( norm(u1-u2) ~= 0 );


%% Problem with breakpoints
dom = [-1 0 pi];
N = chebop(@(x,u) diff(u,2) + cos(x).*u, dom);
N.lbc = 2; 
N.rbc = -1;

x = chebfun(@(x) x, dom);
rhs = sin(x);

%% Try different discretizations
% Start with collocation
pref.discretization = @colloc2;
u3 = solvebvp(N, rhs, pref);

%% Change to ultraS
pref.discretization = @ultraS;
u4 = solvebvp(N, rhs, pref);

%% Did we pass? 
% To pass, both residuals have to be small, but we should not expect u3 and u4
% to be identical!
tol = pref.errTol;
pass(4) = norm(N(u3)-rhs) < tol && ( u3(-1) - 2 < tol) && ( u3(pi) + 1 < tol);
pass(5) = norm(N(u4)-rhs) < tol && ( u4(-1) - 2 < tol) && ( u4(pi) + 1 < tol);
pass(6) = norm(jump(u3,0)) < tol && norm(jump(u4,0)) < tol;
pass(7) = ( norm(u3-u4) ~= 0 );


%% Problem with discontinuous coefficients
dom = [-1 0 pi];
N = chebop(@(x,u) diff(u,2) + abs(x-1).*u, dom);
N.lbc = 2; 
N.rbc = -1;

x = chebfun(@(x) x, dom);
rhs = sin(x);

%% Try different discretizations
% Start with collocation
pref.discretization = @colloc2;
u5 = solvebvp(N, rhs, pref);

%% Change to ultraS
pref.discretization = @ultraS;
u6 = solvebvp(N, rhs, pref);

%% Did we pass? 
% To pass, both residuals have to be small, but we should not expect u3 and u4
% to be identical!
tol = pref.errTol;
pass(8) = norm(N(u5)-rhs) < 10*tol && ...
    ( u5(-1) - 2 < tol) && ( u5(pi) + 1 < tol);
pass(9) = norm(N(u6)-rhs) < 10*tol && ...
    ( u6(-1) - 2 < tol) && ( u6(pi) + 1 < tol);
pass(10) = norm(jump(u5,0)) < tol && norm(jump(u5,0)) < tol;
pass(11) = ( norm(u5-u6) ~= 0 );

end
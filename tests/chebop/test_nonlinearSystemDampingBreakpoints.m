function pass = test_nonlinearSystemDampingBreakpoints(pref)
% Test 2x2 system (sin/cos) where damping is required
%
% Asgeir Birkisson, May 2014.

if ( nargin == 0 )
    pref = cheboppref;
end

tol = 1e-10;

%% Piecewise:
d = [-pi 0 pi];
x = chebfun('x',d);
f = [ 0*x ; 0*x ];

%% COLLOC2

A = chebop(@(x,u,v) [u - diff(v,2); diff(u,2) + cos(v)], d);
A.lbc = @(u,v) [u-1/2; v+1/4];
A.rbc = @(u,v) [v-1/4; u+1/2];

pref.discretization = @colloc2;

[u12, info] = solvebvp(A, f, pref);
u1 = u12{1}; u2 = u12{2};

% Want to check BCs as well.
bcFunLeft = chebfun(A.lbc(u1,u2));
bcFunRight = chebfun(A.rbc(u1,u2));

% And check that we're continuous over breakpoint
u1jump = jump(u1, 0);
u2jump = jump(u2, 0);

pass(1) = abs(info.error) < tol;
pass(2) = norm(bcFunLeft(d(1))) < tol && norm(bcFunRight(d(end))) < tol;
pass(3) = norm(u1jump) < tol && norm(u2jump) < tol;

%% COLLOC1

A = chebop(@(x,u,v) [u - diff(v,2); diff(u,2) + cos(v)], d);
A.lbc = @(u,v) [u-1/2; v+1/4];
A.rbc = @(u,v) [v-1/4; u+1/2];

pref.discretization = @colloc1;

[u34, info] = solvebvp(A, f, pref);
u3 = u34{1}; u4 = u34{2};

% Want to check BCs as well.
bcFunLeft = chebfun(A.lbc(u3,u4));
bcFunRight = chebfun(A.rbc(u3,u4));

% And check that we're continuous over breakpoint
u5jump = jump(u3, 0);
u6jump = jump(u4, 0);

pass(4) = abs(info.error) < tol;
pass(5) = norm(bcFunLeft(d(1))) < tol && norm(bcFunRight(d(end))) < tol;
pass(6) = norm(u5jump) < tol && norm(u6jump) < tol;

%% ULTRAS:

pref.discretization = @ultraS;
[u56, info] = solvebvp(A, f, pref);
u5 = u56{1}; u6 = u56{2};

% Want to check BCs as well.
bcFunLeft = chebfun(A.lbc(u5,u6));
bcFunRight = chebfun(A.rbc(u5,u6));

% And check that we're continuous over breakpoint
u7jump = jump(u5, 0);
u8jump = jump(u6, 0);

pass(7) = abs(info.error) < tol;
pass(8) = norm(bcFunLeft(d(1))) < tol && norm(bcFunRight(d(end))) < tol;
pass(9) = norm(u7jump) < tol && norm(u8jump) < tol;

end
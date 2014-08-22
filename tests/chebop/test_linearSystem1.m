function pass = test_linearSystem1(pref)
% A linear CHEBOP test. This test tests a system of coupled ODEs, both with and
% without breakpoints in the domain of the CHEBOP. Problems are solved using
% chebcolloc1, chebcolloc2 and ultraS discretizations.

% Test 2x2 system (sin/cos)
% Toby Driscoll (extended by AB)

% Note: This test is taken from chebop_systemsolve1 in V4.

if ( nargin == 0 )
    pref = cheboppref;
end

tol = 1e-10;

%% Smooth domain, chebcolloc1:
d = [-pi pi];
A = chebop(@(x,u,v) [u - diff(v) ; diff(u) + v],d);
A.lbc = @(u,v) u + 1;
A.rbc = @(u,v) v;
x = chebfun('x',d);

pref.discretization = @chebcolloc1;
u = mldivide(A, 0, pref);
u1 = u{1}; u2 = u{2};

% Want to check BCs as well.
bcFunLeft = A.lbc(u1,u2);
bcFunRight = A.rbc(u1,u2);

pass(1) = norm( u1 - cos(x),inf) < 100*tol;
pass(2) = norm( u2 - sin(x),inf) < 100*tol;
pass(3) = norm(bcFunLeft(d(1))) < tol && norm(bcFunRight(d(end))) < tol;

%% Smooth domain, chebcolloc2
pref.discretization = @chebcolloc2;
u = mldivide(A, 0, pref);
u1 = u{1}; u2 = u{2};

% Want to check BCs as well.
bcFunLeft = A.lbc(u1,u2);
bcFunRight = A.rbc(u1,u2);

pass(4) = norm( u1 - cos(x),inf) < 100*tol;
pass(5) = norm( u2 - sin(x),inf) < 100*tol;
pass(6) = norm(bcFunLeft(d(1))) < tol && norm(bcFunRight(d(end))) < tol;

%% Smooth domain, ultraS
pref.discretization = @ultraS;
u = mldivide(A, 0, pref);
u1 = u{1}; u2 = u{2};

% Want to check BCs as well.
bcFunLeft = A.lbc(u1,u2);
bcFunRight = A.rbc(u1,u2);

pass(7) = norm( u1 - cos(x),inf) < 100*tol;
pass(8) = norm( u2 - sin(x),inf) < 100*tol;
pass(9) = norm(bcFunLeft(d(1))) < tol && norm(bcFunRight(d(end))) < tol;
%% Piecewise, chebcolloc1:
A.domain = [-pi 0 pi];
pref.discretization = @chebcolloc1;
u = mldivide(A, 0, pref);
u1 = u{1}; u2 = u{2};

% Want to check BCs as well.
bcFunLeft = A.lbc(u1,u2);
bcFunRight = A.rbc(u1,u2);
bcValLeft = bcFunLeft(d(1));
bcValRight = bcFunRight(d(2));

% And check that we're continuous over breakpoint
u1jump = jump(u1, 0);
u2jump = jump(u2, 0);

pass(10) = norm( u1 - cos(x),inf) < 2000*tol;
pass(11) = norm( u2 - sin(x),inf) < 2000*tol;
pass(12) = norm(bcValLeft) < tol && norm(bcValRight) < tol;
pass(13) = norm(u1jump) < tol && norm(u2jump) < tol;

%% Piecewise, chebcolloc2:
pref.discretization = @chebcolloc2;
u = mldivide(A, 0, pref);
u1 = u{1}; u2 = u{2};

% Want to check BCs as well.
bcFunLeft = A.lbc(u1,u2);
bcFunRight = A.rbc(u1,u2);
bcValLeft = bcFunLeft(d(1));
bcValRight = bcFunRight(d(2));

% And check that we're continuous over breakpoint
u1jump = jump(u1, 0);
u2jump = jump(u2, 0);

pass(14) = norm( u1 - cos(x),inf) < 2000*tol;
pass(15) = norm( u2 - sin(x),inf) < 2000*tol;
pass(16) = norm(bcValLeft) < tol && norm(bcValRight) < tol;
pass(17) = norm(u1jump) < tol && norm(u2jump) < tol;

%% Piecewise, ultraS:
pref.discretization = @ultraS;
u = mldivide(A, 0, pref);
u1 = u{1}; u2 = u{2};

% Want to check BCs as well.
bcFunLeft = A.lbc(u1,u2);
bcFunRight = A.rbc(u1,u2);
bcValLeft = bcFunLeft(d(1));
bcValRight = bcFunRight(d(2));

% And check that we're continuous over breakpoint
u1jump = jump(u1, 0);
u2jump = jump(u2, 0);

pass(18) = norm( u1 - cos(x),inf) < 2000*tol;
pass(19) = norm( u2 - sin(x),inf) < 2000*tol;
pass(20) = norm(bcValLeft) < tol && norm(bcValRight) < tol;
pass(21) = norm(u1jump) < tol && norm(u2jump) < tol;
end



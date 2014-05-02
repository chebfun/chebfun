function pass = test_linearSystem1(pref)
% Test 2x2 system (sin/cos)
% Toby Driscoll (extended by AB)

% Note: This test is taken from chebop_systemsolve1 in V4.

if ( nargin == 0 )
    pref = cheboppref;
end

tol = 1e-10;

pref.display = 'iter';

%% Smooth domain:
d = [-pi pi];
A = chebop(@(x,u,v) [u - diff(v) ; diff(u) + v],d);
A.lbc = @(u,v) u + 1;
A.rbc = @(u,v) v;
x = chebfun('x',d);
f = [ 0*x ; 0*x ];
u = mldivide(A, f, pref);

u1 = u{1}; u2 = u{2};

% Want to check BCs as well.
bcFunLeft = A.lbc(u1,u2);
bcFunRight = A.rbc(u1,u2);

pass(1) = norm( u1 - cos(x),inf) < 100*tol;
pass(2) = norm( u2 - sin(x),inf) < 100*tol;
pass(3) = norm(bcFunLeft(d(1))) < tol && norm(bcFunRight(d(end))) < tol;

%% Piecewise:
A.domain = [-pi 0 pi];
u = mldivide(A, f, pref);
u1 = u{1}; u2 = u{2};

% Want to check BCs as well.
bcFunLeft = A.lbc(u1,u2);
bcFunRight = A.rbc(u1,u2);
bcValLeft = bcFunLeft(d(1));
bcValRight = bcFunRight(d(2));

% And check that we're continuous over breakpoint
u1jump = jump(u1, 0);
u2jump = jump(u2, 0);

pass(4) = norm( u1 - cos(x),inf) < 2000*tol;
pass(5) = norm( u2 - sin(x),inf) < 2000*tol;
pass(6) = norm(bcValLeft) < tol && norm(bcValRight) < tol;
pass(7) = norm(u1jump) < tol && norm(u2jump) < tol;

end



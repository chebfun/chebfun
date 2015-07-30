function pass = test_nonlinSys1Breaks_US(pref)
% Test 2x2 system (sin/cos). This is pecewiseificaion of the test
%       test_nonlinearSystem1
%
% Asgeir Birkisson, April 2014.

if ( nargin == 0 )
    pref = cheboppref;
end

tol = 1e-10;

d = [-pi 0 pi];
x = chebfun('x',d);
f = [ 0*x ; 0*x ];

%% Piecewise (ultraS):
pref.discretization = @ultraS;

A = chebop(@(x,u,v) [u - diff(v,2) + u.^2; diff(u) + sin(v)],d);
A.lbc = @(u,v) u-1;
A.rbc = @(u,v) [v-1/2; diff(v)];

u = mldivide(A, f, pref);
u1 = u{1}; u2 = u{2};

% Want to check BCs as well.
bcFunLeft = A.lbc(u1,u2);
bcFunRight = chebfun(A.rbc(u1,u2));

% And check that we're continuous over breakpoint
u1jump = jump(u1, 0);
u2jump = jump(u2, 0);

pass(1) = norm( chebfun(A(x, u1, u2))) < tol;
pass(2) = norm(bcFunLeft(d(1))) < tol && norm(bcFunRight(d(end))) < tol;
pass(3) = norm(u1jump) < tol && norm(u2jump) < tol;

end

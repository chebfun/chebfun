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

A = chebop(@(x,u,v) [u - diff(v,2) + (1-x.^2).*u.^2; diff(u,2) + sin(v)], d);
A.lbc = @(u,v) [u-1; v+1];
A.rbc = @(u,v) [v-1/2; diff(v)];

pref.discretization = @colloc2;
u56 = mldivide(A, f, pref);
u5 = u56{1}; u6 = u56{2};

% Want to check BCs as well.
bcFunLeft = chebfun(A.lbc(u5,u6));
bcFunRight = chebfun(A.rbc(u5,u6));

% And check that we're continuous over breakpoint
u5jump = jump(u5, 0);
u6jump = jump(u6, 0);

pass(1) = norm( chebfun(A(x, u5, u6))) < tol;
pass(2) = norm(bcFunLeft(d(1))) < tol && norm(bcFunRight(d(end))) < tol;
pass(3) = norm(u5jump) < tol && norm(u6jump) < tol;

%% Piecewise:

pref.discretization = @ultraS;
u78 = mldivide(A, f, pref);
u7 = u78{1}; u8 = u78{2};

% Want to check BCs as well.
bcFunLeft = chebfun(A.lbc(u7,u8));
bcFunRight = chebfun(A.rbc(u7,u8));

% And check that we're continuous over breakpoint
u7jump = jump(u7, 0);
u8jump = jump(u8, 0);

pass(4) = norm( chebfun(A(x, u7, u8))) < tol;
pass(5) = norm(bcFunLeft(d(1))) < tol && norm(bcFunRight(d(end))) < tol;
pass(6) = norm(u7jump) < tol && norm(u8jump) < tol;

end
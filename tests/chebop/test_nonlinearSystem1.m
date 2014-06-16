function pass = test_nonlinearSystem1(pref)
% Test 2x2 system (sin/cos). This is nonlinearification of the test
%       test_linearSystem1
%
% Asgeir Birkisson, April 2014.

if ( nargin == 0 )
    pref = cheboppref;
end

tol = 1e-10;

% Smooth domain:
d = [-pi pi];
x = chebfun('x',d);
f = [ 0*x ; 0*x ];

%% Colloc2
pref.discretization = @colloc2;

A = chebop(@(x,u,v) [u - diff(v,2) + u.^2; diff(u) + sin(v)],d);
A.lbc = @(u,v) u-1;
A.rbc = @(u,v) [v-1/2; diff(v)];

u12 = mldivide(A, f, pref);
u1 = u12{1}; u2 = u12{2};

% Want to check BCs as well.
bcFunLeft = A.lbc(u1,u2);
bcFunRight = chebfun(A.rbc(u1,u2));

pass(1) = norm( chebfun(A(x, u1, u2))) < tol;
pass(2) = norm(bcFunLeft(d(1))) < tol && norm(bcFunRight(d(end))) < tol;

%% Try with ultraS
pref.discretization = @ultraS;
u34 = mldivide(A, f, pref);

u3 = u34{1}; u4 = u34{2};

% Want to check BCs as well.
bcFunLeft = A.lbc(u3,u4);
bcFunRight = chebfun(A.rbc(u3,u4));

pass(3) = norm( chebfun(A(x, u3, u4))) < tol;
pass(4) = norm(bcFunLeft(d(1))) < tol && norm(bcFunRight(d(end))) < tol;

%% Try with colloc1
pref.discretization = @colloc1;
u56 = mldivide(A, f, pref);

u5 = u56{1}; u6 = u56{2};

% Want to check BCs as well.
bcFunLeft = A.lbc(u5,u6);
bcFunRight = chebfun(A.rbc(u5,u6));

pass(5) = norm( chebfun(A(x, u5, u6))) < tol;
pass(6) = norm(bcFunLeft(d(1))) < tol && norm(bcFunRight(d(end))) < tol;

%% We expect the solutions to be similar, but not identical!
pass(7) = ( (norm(u12 - u34) ~= 0) && (norm(u12 - u56) ~= 0) && ...
    (norm(u34 - u56) ~= 0));

end



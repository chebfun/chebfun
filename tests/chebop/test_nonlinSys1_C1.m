function pass = test_nonlinSys1_C1(pref)
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

%% Colloc1
pref.discretization = @chebcolloc1;

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

end



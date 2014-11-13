function pass = test_nonlinSys2_C1(pref)
% Test 2x2 system (sin/cos). This is the same problem as test_nonlinSystem1,
% but using CHEBMATRIX {} syntax.
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

A = chebop(@(x,u) [u{1} - diff(u{2}, 2) + u{1}.^2; diff(u{1}) + sin(u{2})], d);
A.lbc = @(u) u{1} - 1;
A.rbc = @(u) [u{2} - 1/2; diff(u{2})];

u = mldivide(A, f, pref);

% Check the residual
err1 = norm(A(u) - f );

% Want to check BCs as well.
bcFunLeft = A.lbc(u);
bcFunRight = chebfun(A.rbc(u));
err2 = [norm(bcFunLeft(d(1))), norm(bcFunRight(d(end)))];

pass(1) = err1 < tol;
pass(2) = all( err2 < tol );

end



function pass = test_nonlinearSystem2(pref)
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

%% Colloc2
pref.discretization = @colloc2;

A = chebop(@(x,u) [u{1} - diff(u{2}, 2) + u{1}.^2; diff(u{1}) + sin(u{2})], d);
A.lbc = @(u) u{1} - 1;
A.rbc = @(u) [u{2} - 1/2; diff(u{2})];

u12 = mldivide(A, f, pref);

% Check the residual
err1 = norm(A(u12) - f );

% Want to check BCs as well.
bcFunLeft = A.lbc(u12);
bcFunRight = chebfun(A.rbc(u12));
err2 = [norm(bcFunLeft(d(1))), norm(bcFunRight(d(end)))];

pass(1) = err1 < tol;
pass(2) = all( err2 < tol );
%% Try with ultraS
pref.discretization = @ultraS;
u34 = mldivide(A, f, pref);

% Check the residual
err3 = norm(A(u34) - f );

% Want to check BCs as well.
bcFunLeft = A.lbc(u34);
bcFunRight = chebfun(A.rbc(u34));
err4 = [norm(bcFunLeft(d(1))), norm(bcFunRight(d(end)))];

pass(3) = err3 < tol;
pass(4) = all( err4 < tol );
%% Try with colloc1
pref.discretization = @colloc1;
u56 = mldivide(A, f, pref);

% Check the residual
err5 = norm(A(u56) - f );

% Want to check BCs as well.
bcFunLeft = A.lbc(u56);
bcFunRight = chebfun(A.rbc(u56));
err6 = [norm(bcFunLeft(d(1))), norm(bcFunRight(d(end)))];

pass(5) = err5 < tol;
pass(6) = all( err6 < tol );

%% We expect the solutions to be similar, but not identical!
pass(7) = ( (norm(u12 - u34) ~= 0) && (norm(u12 - u56) ~= 0) && ...
    (norm(u34 - u56) ~= 0));

end



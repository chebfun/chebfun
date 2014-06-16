function pass = test_linearSystem2(pref)
% Test solution of a 2x2 system
% Toby Driscoll (extended by AB)

% Note: This test is taken from chebop_systemsolve1 in V4.

%% Setup:
if ( nargin == 0 )
    pref = cheboppref;
end

tol = 1e-10;

%% Smooth domain:
d = [-1 1];
A = chebop(@(x, u) [diff(u{1}) + u{1} + 2*u{2} ; ...
    diff(u{1}) - u{1} + diff(u{2})], d);
A.lbc = @(u) u{1}+diff(u{1});
A.rbc = @(u) diff(u{2});
x = chebfun('x',d);

f = [ exp(x) ; chebfun(1,d) ];
u = mldivide(A, f, pref);

% Residual of differential equation:
err1 = norm(A(u) - f);

% Want to check BCs as well.
bcFunLeft = A.lbc(u);
bcFunRight = A.rbc(u);
err2 = [norm(bcFunLeft(d(1))), norm(bcFunRight(d(end)))];

% Happy?
pass(1) = norm( err1 ) < tol;
pass(2) = all(err2 < tol);
%% Piecewise:
A.domain = [-1 0 1];
u = mldivide(A, f, pref);

% Norm of the residual:
err3 = norm( A(u) - f );

% Want to check BCs as well:
bcFunLeft = A.lbc(u);
bcFunRight = A.rbc(u);
err4 = [norm(bcFunLeft(d(1))), norm(bcFunRight(d(end)))];

% And check that we're continuous over the breakpoint:
u1jump = jump(u{1}, 0);
u2jump = jump(u{2}, 0);
err5 = [norm(u1jump), norm(u2jump)];

% Happy?
pass(3) = err3 < tol;
pass(4) = all(err4 < tol);
pass(5) = all(err5 < tol);

%% Now use a function handle, rathen than anonymous function.

% Let's change the discretization as well, just because we can:
pref.discretization = @ultraS;
A = chebop(@myop, d);
A.lbc = @(u) u{1}+diff(u{1});
A.rbc = @(u) diff(u{2});

% Here we need to pass information about the number of variables that A
% operates on:
A.numVars = 2;
u = mldivide(A, f, pref);

% Residual of differential equation:
err6 = norm(A(u) - f);

% Want to check BCs as well.
bcFunLeft = A.lbc(u);
bcFunRight = A.rbc(u);
err7 = [norm(bcFunLeft(d(1))), norm(bcFunRight(d(end)))];

% Happy?
pass(6) = norm( err6 ) < tol;
pass(7) = all(err7 < tol);

%% Now try piecewise again!
A.domain = [-1 0 1];
u = mldivide(A, f, pref);

% Norm of the residual:
err8 = norm( A(u) - f );

% Want to check BCs as well:
bcFunLeft = A.lbc(u);
bcFunRight = A.rbc(u);
err9 = [norm(bcFunLeft(d(1))), norm(bcFunRight(d(end)))];

% And check that we're continuous over the breakpoint:
u1jump = jump(u{1}, 0);
u2jump = jump(u{2}, 0);
err10 = [norm(u1jump), norm(u2jump)];

% Happy?
pass(8)  = err8 < tol;
pass(9)  = all(err9 < tol);
pass(10) = all(err10 < tol);

%% Finally, let's include x in the anonmoys function argument list:
d = [-1 1];
A = chebop(@(u) [diff(u{1}) + u{1} + 2*u{2} ; ...
    diff(u{1}) - u{1} + diff(u{2})], d);
A.lbc = @(u) u{1}+diff(u{1});
A.rbc = @(u) diff(u{2});

% And for fun, let's try another discretization
pref.discretization = @ultraS;
u = mldivide(A, f, pref);

% Residual of differential equation:
err11 = norm(A(u) - f);

% Want to check BCs as well.
bcFunLeft = A.lbc(u);
bcFunRight = A.rbc(u);
err12 = [norm(bcFunLeft(d(1))), norm(bcFunRight(d(end)))];

% Happy?
pass(11) = norm( err11 ) < tol;
pass(12) = all(err12 < tol);

end

function out = myop(x, u)

out = [diff(u{1}) + u{1} + 2*u{2} ; ...
    diff(u{1}) - u{1} + diff(u{2})];

end

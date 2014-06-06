function pass = test_linearSystem2(pref)
% Test solution of a 2x2 system
% Toby Driscoll (extended by AB)

% Note: This test is taken from chebop_systemsolve1 in V4.

if ( nargin == 0 )
    pref = cheboppref;
end
%%
tol = 1e-10;

% Smooth domain:
d = [-1 1];
A = chebop(@(x, u) [diff(u{1}) + u{1} + 2*u{2} ; ...
    diff(u{1}) - u{1} + diff(u{2})], d);
A.lbc = @(u) u{1}+diff(u{1});
A.rbc = @(u) diff(u{2});
x = chebfun('x',d);

% TODO: Shouldn't need to pass an initial guess, but necessary to get the
% correct dimension information!
A.init = [0*x; 0*x];
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
end

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
f = [ exp(x) ; chebfun(1,d) ];
u = mldivide(A, f, pref);

u1 = u{1}; u2 = u{2};
pass(1) = norm( diff(u1)+u1+2*u2-exp(x)) < tol;
pass(2) = norm( diff(u1)-u1+diff(u2)-1 ) < tol;

% Want to check BCs as well.
bcFunLeft = A.lbc(u1,u2);
bcFunRight = A.rbc(u1,u2);
pass(3) = norm(bcFunLeft(d(1))) < tol && norm(bcFunRight(d(end))) < tol;

%% Piecewise:
A.domain = [-1 0 1];
u = mldivide(A, f, pref);
u1 = u{1}; u2 = u{2};

err1 = diff(u1) + u1 + 2*u2 - exp(x);
err2 = diff(u1) - u1 + diff(u2) - 1;

% Want to check BCs as well.
bcFunLeft = A.lbc(u1,u2);
bcFunRight = A.rbc(u1,u2);
bcValLeft = bcFunLeft(d(1));
bcValRight = bcFunRight(d(2));

% And check that we're continuous over breakpoint
u1jump = jump(u1, 0);
u2jump = jump(u2, 0);

pass(4) = norm( err1 ) < tol;
pass(5) = norm( err2 ) < tol;
pass(6) = norm(bcValLeft) < tol && norm(bcValRight) < tol;
pass(7) = norm(u1jump) < tol && norm(u2jump) < tol;

end

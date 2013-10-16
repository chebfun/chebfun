% Test QR factorization of Chebfun quasimatrices.

function pass = test_qr(pref)

% [TODO]: This test needs to be much more extensive.

if ( nargin == 0 )
    pref = chebfun.pref();
end

%% Test a smooth CHEBFUN:
A = chebfun(@(x) [x 1i*x 1+0*x 1+1i+0*x (2-1i)*x], [0, 1], pref);
pass(1) = (rank(A) == 2);
[Q, R] = qr(A);

tol = epslevel(A);
% Test 2 confirms that the columns of Q are orthonormal:
pass(2) = abs(cond(Q) - 1) < 1e-13*(tol/eps);

% Test 3 confirms that the Q and R factors have the right product.
pass(3) = norm(A - Q*R) < 1e-13*(tol/eps);

%% Test with a breakpoint:
A(.5,:) = feval(A, .5);
pass(4) = (rank(A) == 2);
[Q, R] = qr(A);

tol = epslevel(A);
% Test 2 confirms that the columns of Q are orthonormal:
pass(5) = abs(cond(Q) - 1) < 1e-13*(tol/eps);

% Test 3 confirms that the Q and R factors have the right product.
pass(6) = norm(A - Q*R) < 1e-13*(tol/eps);

% Check for bug in piecewise QR from commit 9ba78c2a.
f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 0 1], pref);
[Q, R] = qr(f);
pass(7) = norm(f - Q*R) < 10*vscale(f)*epslevel(f);

% Check QR of a CHEBFUN with one column and breakpoints.
f = chebfun(@(x) 1 + 0*x, [-1 0 1], pref);
[Q, R] = qr(f);
pass(8) = norm(f - Q*R) < 10*vscale(f)*epslevel(f) && ...
    abs(Q'*Q - 1) < 10*vscale(Q)*epslevel(Q);

% Check QR of a CHEBFUN based on CHEBTECH1.
p = pref;
p.chebfun.tech = 'chebtech1';
f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 0 1], p);
[Q, R] = qr(f);
pass(9) = strcmp(class(f.funs{1}.onefun), 'chebtech1') && ...
    norm(f - Q*R) < 10*vscale(f)*epslevel(f) && ...
    norm(Q'*Q - eye(3), 'fro') < 10*vscale(Q)*epslevel(Q);

end

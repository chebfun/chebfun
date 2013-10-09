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

end
% Test for bary.m.

function pass = test_bary(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

xr_row = xr.';
xr_mtx = reshape(xr, [100 10]);
xr_3mtx = reshape(xr, [10 10 10]);

% Set an error tolerance.
tol = 1.0e-14;

% Check evaluation for a scalar polynomial.
p = @(x) x.^4 - 2*x.^3 + 3*x.^2 - 2*x + 1;

[xk, ignored, vk] = chebpts(16);
err = norm(p(xr) - bary(xr, p(xk), xk, vk), Inf);
pass(1) = err < tol;

err = norm(p(xr(1)) - bary(xr(1), p(xk), xk, vk), Inf);
pass(2) = err < tol;

[xk, ignored, vk] = chebpts(16, 1);
err = norm(p(xr) - bary(xr, p(xk), xk, vk), Inf);
pass(3) = err < tol;

err = norm(p(xr(1)) - bary(xr(1), p(xk), xk, vk), Inf);
pass(4) = err < tol;

% Check row vector and matrix input with scalar polynomial.
[xk, ignored, vk] = chebpts(16);

difference = p(xr_row) - bary(xr_row, p(xk), xk, vk);
err = norm(difference(:), Inf);
pass(5) = err < tol;

difference = p(xr_mtx) - bary(xr_mtx, p(xk), xk, vk);
err = norm(difference(:), Inf);
pass(6) = err < tol;

difference = p(xr_3mtx) - bary(xr_3mtx, p(xk), xk, vk);
err = norm(difference(:), Inf);
pass(7) = err < tol;

% Check evaluation for an array of two polynomials.
q1 = @(x) x.^4 - 2*x.^3 + 3*x.^2 - 2*x + 1;
q2 = @(x) 2*x.^4 - 3*x.^3 + x.^2 - 5*x + 2;
q = @(x) [q1(x) q2(x)];

[xk, ignored, vk] = chebpts(16);
difference = q(xr) - bary(xr, q(xk), xk, vk);
err = norm(difference(:), Inf);
pass(8) = err < tol;

difference = q(xr(1)) - bary(xr(1), q(xk), xk, vk);
err = norm(difference(:), Inf);
pass(9) = err < tol;

[xk, ignored, vk] = chebpts(16, 1);
difference = q(xr) - bary(xr, q(xk), xk, vk);
err = norm(difference(:), Inf);
pass(10) = err < tol;

difference = q(xr(1)) - bary(xr(1), q(xk), xk, vk);
err = norm(difference(:), Inf);
pass(11) = err < tol;

% Check row vector and matrix input with array of polynomials.
[xk, ignored, vk] = chebpts(16);

difference = q(xr_row) - bary(xr_row, q(xk), xk, vk);
err = norm(difference(:), Inf);
pass(12) = err < tol;

difference = q(xr_mtx) - bary(xr_mtx, q(xk), xk, vk);
err = norm(difference(:), Inf);
pass(13) = err < tol;

% Make sure (a smaller version of) the example from the file doesn't crash.
try
    x = chebpts(181);
    f = 1./(1 + 25*x.^2);
    xx = linspace(-1, 1, 100);
    [xx, yy] = meshgrid(xx, xx);
    ff = bary(xx + 1i*yy, f);
    pass(14) = true;
catch ME
    pass(14) = false;
end

end

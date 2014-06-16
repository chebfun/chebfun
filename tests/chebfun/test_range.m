% Test file for @chebfun/range.m.

function pass = test_range(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Empty input.
pass(1) = isempty(range(chebfun()));

% Scalar row and column input.
f = chebfun(@(x) x, [-1 0 1]);
err = abs(range(f) - 2);
tol = 10*vscale(f)*epslevel(f);
pass(2) = err < tol;

ft = f.';
err = abs(range(ft) - 2);
tol = 10*vscale(ft)*epslevel(ft);
pass(3) = err < tol;

% Array-valued column and column quasimatrix inputs.
fa = chebfun(@(x) [x x.^2], [0 1 2]);
fq = cheb2quasi(fa);

err = norm(range(fa) - [2 4], Inf);
tol = 10*vscale(fa)*epslevel(fa);
pass(4) = err < tol;

err = norm(range(fq) - [2 4], Inf);
tol = 10*vscale(fq)*epslevel(fq);
pass(5) = err < tol;

ra = range(fa, 2);
err = norm(feval(ra, x) - r_exact(x), Inf);
tol = 10*vscale(ra)*epslevel(ra);
pass(6) = err < tol;

rq = range(fq, 2);
err = norm(feval(rq, x) - r_exact(x), Inf);
tol = 10*vscale(rq)*epslevel(rq);
pass(7) = err < tol;

% Array-valued row and row quasimatrix inputs.
fat = fa.';
fqt = fq.';

err = norm(range(fat, 2) - [2 ; 4], Inf);
tol = 10*vscale(fa)*epslevel(fat);
pass(8) = err < tol;

err = norm(range(fqt, 2) - [2 ; 4], Inf);
tol = 10*vscale(fq)*epslevel(fqt);
pass(9) = err < tol;

rat = range(fat, 1);
err = norm(feval(rat, x) - r_exact(x).', Inf);
tol = 10*vscale(rat)*epslevel(rat);
pass(10) = err < tol;

% Dimension >= 3 input.
pass(11) = iszero(range(f, 3));

end

function y = r_exact(x)
    y = zeros(size(x));
    ind1 = x <= 1;
    y(ind1) = x(ind1) - x(ind1).^2;
    ind2 = x > 1;
    y(ind2) = x(ind2).^2 - x(ind2);
end

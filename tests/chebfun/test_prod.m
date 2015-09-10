% Test file for chebfun/prod.m.

function pass = test_prod(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

tol = 100*eps;

f = chebfun(@(x) x + 2, pref);
fa = [f f];

% Scalar-valued.
err = abs(prod(f) - 27*exp(-2));
pass(1) = err < tol;
err = abs(prod(f.') - 27*exp(-2));
pass(2) = err < tol;

% Array-valued, continuous dimension.
err = norm(abs(prod(fa) - [27*exp(-2) 27*exp(-2)]), Inf);
pass(3) = err < tol;
pass(4) = isequal(prod(fa), prod(fa, 1));

% Array-valued, continuous dimension, transposed.
err = norm(abs(prod(fa.') - [27*exp(-2) 27*exp(-2)].'), Inf);
pass(5) = err < tol;
pass(6) = isequal(prod(fa.'), prod(fa.', 2));

fq = cheb2quasi(fa);

% Quasimatrix, continuous dimension.
err = norm(abs(prod(fq) - [27*exp(-2) 27*exp(-2)]), Inf);
pass(7) = err < tol;
pass(8) = isequal(prod(fq), prod(fq, 1));

% Quasimatrix, continuous dimension, transposed.
err = norm(abs(prod(fq.') - [27*exp(-2) 27*exp(-2)].'), Inf);
pass(9) = err < tol;
pass(10) = isequal(prod(fq.'), prod(fq.', 2));

f = chebfun(@(x) x, pref);
fa = [f f];

% Array-valued, discrete dimension.
pfa2 = prod(fa, 2);
pass(11) = ~pfa2.isTransposed;
err = norm(f.^2 - pfa2, Inf);
pass(12) = err < tol;

% Array-valued, discrete dimension, transposed.
pfat2 = prod(fa.', 1);
pass(13) = pfat2.isTransposed;
err = norm((f.').^2 - pfat2, Inf);
pass(14) = err < tol;

fq = cheb2quasi(fa);

% Quasimatrix, discrete dimension.
pfq2 = prod(fq, 2);
pass(15) = ~pfq2.isTransposed;
err = norm(f.^2 - pfq2, Inf);
pass(16) = err < tol;

% Quasimatrix, discrete dimension, transposed.
pfqt2 = prod(fq.', 1);
pass(17) = pfat2.isTransposed;
err = norm((f.').^2 - pfqt2, Inf);
pass(18) = err < tol;

end

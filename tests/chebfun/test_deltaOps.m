% Test some operations involving chebfuns with delta functions.

function pass = test_deltaOps(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

tol = pref.deltaPrefs.deltaTol;

% These tests used to live in the scripts in tests/deltafun, but since they
% create chebfuns, they are better placed here.

% Some tests for cumsum() based on examples provided by LNT.
n = 6;
x = chebfun('x',[0 n]);
f = 0.5*sin(x);
for j = 1:n-1
  f = f + randn*dirac(x-j);
end
err = norm(diff(cumsum(f)) - f);
pass(1) = err < tol;

x = chebfun('x');
f = dirac(x-.5) + dirac(x) + dirac(x+.5) + heaviside(x);
err = norm(diff(cumsum(f)) - f);
pass(2) = err < tol;

x = chebfun('x');
f = sign(x)+sign(x-.5);
f2 = f(-1) + cumsum(diff(f));
err = norm(f - f2);
pass(3) = err < tol;

% A test for diff() based on an example from LNT.
savedPrefs = chebfunpref();
chebfunpref.setDefaults('enableDeltaFunctions', true);
try
    x = chebfun('x', [0 5]);
    f = 0.5*sin(x);
    A = randn(4, 1);
    for j = 1:4
      f = f + A(j)*dirac(x-j);
    end
    F = cumsum(.5*sin(x));
    for j = 1:4
      F = F + A(j)*heaviside(x-j);
    end
    err = norm(diff(F) - (f - f(0)));
    pass(4) = err < tol;
catch ME
    chebfunpref.setDefaults(savedPrefs);
    rethrow(ME)
end

chebfunpref.setDefaults(savedPrefs);

% A test involving innerProduct:
x = chebfun('x');
pass(5) = innerProduct(dirac(x), x) < eps;

% Test quasi matrix construction
x = chebfun('x');
A = [];
A = [ A, dirac(x-1) ];
A = [ A, dirac(x-0.5) ];
A = [ A, dirac(x-0.0) ];
pass(6) = (norm(A(:, 2) - dirac(x-.5)) < eps) && (size(A,2) == 3);

end

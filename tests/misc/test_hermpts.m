function pass = test_hermpts(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% TODO: These values were computed in V4. Should perhaps check more carefully?

% Choose a tolerance:
tol = 10*pref.eps;

% Test a small n (using REC)
n = 42;
[x] = hermpts(n);
pass(1) = all(size(x) == [n, 1]);
[x, w, v] = hermpts(n);
pass(2) = all(size(x) == [n, 1]) && all(size(w) == [1, n]) && ...
    all(size(v) == [n, 1]);
pass(3) = abs(w*x) < tol && abs(w*x.^2 - sqrt(pi)/2) < tol;
pass(4) = abs(x(37) - 5.660357581283058) < 10*tol;
pass(5) = abs(w(17) - 0.032202101288908) < tol;
pass(6) = abs(v(17) - 0.311886101735772) < tol;

% Test a larger n (using ASY)
n = 251;
[x] = hermpts(n);
pass(11) = all(size(x) == [n, 1]);
[x, w, v] = hermpts(n);
pass(7) = all(size(x) == [n, 1]) && all(size(w) == [1, n]) && ...
    all(size(v) == [n, 1]);
pass(8) = abs(w*x) < tol && abs(w*x.^2 - sqrt(pi)/2) < 300*tol;
pass(9) = abs(x(37) - -13.292221459334638) < 4*tol;
pass(10) = abs(w(123) - 0.117419270715955) < 10*tol;
pass(11) = abs(v(123) - 0.915560323259764) < 100*tol;

% Test on 'prob':
[x2, w2, v2] = hermpts(n, 'prob');
pass(12) = norm(x - x2./sqrt(2), inf) < tol;
pass(13) = norm(w - w2./sqrt(2), inf) < tol;

[x3, w3, v3] = hermpts(n, 'phys');
pass(14) = norm(x - x3, inf) == 0;
pass(15) = norm(w - w3, inf) == 0;

end
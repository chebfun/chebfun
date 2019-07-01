function pass = test_lagpts(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% TODO: These values were computed in V4. Should perhaps check more carefully?

% Choose a tolerance:
tol = 1e3*pref.chebfuneps;

% Test a small n (using GQ)
n = 42;
[x] = lagpts(n);
pass(1) = all(size(x) == [n, 1]);
[x, w, v] = lagpts(n);
pass(2) = all(size(x) == [n, 1]) && all(size(w) == [1, n]) && ...
    all(size(v) == [n, 1]);

pass(3) = (abs(w*x - 1) <= tol) && (abs(w*x.^2 - 2) <= tol);
pass(4) = abs(x(37) - 98.388267163326702) < tol;
pass(5) = abs(w(7) - 0.055372813167092) < tol;
pass(6) = abs(v(17) - 0.002937421407003) < tol;

% Test a larger n (using ASY)
n = 251;
[x] = lagpts(n);
pass(7) = all(size(x) == [n, 1]);
[x, w, v] = lagpts(n);
pass(8) = all(size(x) == [n, 1]) && all(size(w) == [1, n]) && ...
    all(size(v) == [n, 1]);
pass(9) = abs(w*x - 1) < 100*tol && abs(w*x.^2 - 2) < 100*tol;
pass(10) = abs(x(37) - 13.309000189442097) < tol;
pass(11) = abs(w(3) - 0.050091759039996) < tol;
pass(12) = abs(v(3) - 0.214530194346947) < 10*tol;

% Test a different interval (using GQ)
n = 42;
[x, w, v] = lagpts(n, [1, inf]);
e = exp(1);
pass(13) = abs(w*x - 2/e) < tol && abs(w*x.^2 - 5/e) < tol;

% Put the inf on the left:
[x, w, v] = lagpts(n, [-inf, -1]);
e = exp(1);
pass(14) = abs(w*x + 2/e) < tol && abs(w*x.^2 - 5/e) < tol;

end


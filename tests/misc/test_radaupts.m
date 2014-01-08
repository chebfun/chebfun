function pass = test_lobpts(pref)

% Choose a tolerance:
tol = 1e-14;

n = 42;
[x] = radaupts(n);
pass(1) = all(size(x) == [n, 1]);
[x, w, v] = radaupts(n);
pass(2) = all(size(x) == [n, 1]) && all(size(w) == [1, n]) && all(size(v) == [n, 1]);
pass(3) = abs(w*x) < tol && abs(w*x.^2 - 2/3) < tol;
pass(4) = abs(x(37) - 0.908847278001044) < tol;
pass(5) = abs(w(37) - 0.031190846817016) < tol;
pass(6) = abs(v(37) + 0.171069152683909) < tol;
pass(7) = x(1) == -1;

end
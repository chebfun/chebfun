function pass = test_legpts(pref)

% Choose a tolerance:
tol = 1e-14;

% Test a small n (using REC)
n = 42;
[x] = legpts(n);
pass(1) = all(size(x) == [n, 1]);
[x, w, v] = legpts(n);
pass(2) = all(size(x) == [n, 1]) && all(size(w) == [1, n]) && all(size(v) == [n, 1]);
pass(4) = abs(w*x) < tol && abs(w*x.^2 - 2/3) < tol;
pass(4) = abs(x(37) - 0.910959724904127) < tol;
pass(5) = abs(w(37) - 0.030479240699603) < tol;
pass(6) = abs(v(37) - 0.265155501739424) < tol;

% Test on [0, 10]:
[x, w, v] = legpts(n, [0, 10]);
pass(7) = abs(w*x - 50) < 10*tol && abs(w*x.^2 - 1000/3) < 100*tol;
pass(8) = abs(x(38) - 9.694617786774941) < tol;
pass(9) = abs(w(38) - 0.127114797630565) < tol;
pass(10) = abs(v(38) + 0.202027188941007) < tol;

% Test a larger n (using ASY)
n = 251;
[x] = legpts(n);
pass(11) = all(size(x) == [n, 1]);
[x, w, v] = legpts(n);
pass(12) = all(size(x) == [n, 1]) && all(size(w) == [1, n]) && all(size(v) == [n, 1]);
pass(13) = abs(w*x) < tol && abs(w*x.^2 - 2/3) < tol;
pass(14) = abs(x(37) + 0.896467746955729) < tol;
pass(15) = abs(w(37) - 0.005535005742012) < tol;
pass(16) = abs(v(37) - 0.294960654628873) < tol;

% Test on [0, 10]:
[x, w, v] = legpts(n, [0, 10]);
pass(17) = abs(w*x - 50) < 10*tol && abs(w*x.^2 - 1000/3) < 100*tol;
pass(18) = abs(x(38) - 0.545685271938239) < tol;
pass(19) = abs(w(38) - 0.028372255931272) < tol;
pass(20) = abs(v(38) + 0.306176997099458) < tol;

end
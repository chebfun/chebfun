function pass = test_jacpts(pref)

% Choose a tolerance:
tol = 1e-14;

% Test a small n (using REC)
n = 42;
a = -.1; b = .3;
x = jacpts(n, a, b);
pass(1) = all(size(x) == [n, 1]);
[x, w, v] = jacpts(n, a, b);
pass(2) = all(size(x) == [n, 1]) && all(size(w) == [1, n]) && all(size(v) == [n, 1]);
pass(3) = abs(w*x - 0.363593965943934) < tol && abs(w*x.^2 - 0.670376374709129) < tol;
pass(4) = abs(x(37) - 0.912883347814032) < tol;
pass(5) = abs(w(37) - 0.046661910947553) < tol;
pass(6) = abs(v(37) - 0.320696510075909) < tol;

% Test on [0, 10]:
[x, w, v] = jacpts(n, a, b, [0, 10]);
pass(7) = abs(w*x - 81.519974175439359) < 10*tol && abs(w*x.^2 - 585.9248143859702) < 100*tol;
pass(8) = abs(x(38) - 9.702339316456870) < tol;
pass(9) = abs(w(38) - 0.279566831611687) < tol;
pass(10) = abs(v(38) + 0.248832970298009) < tol;

% Test a larger n (using ASY)
a = -.7; b = 1.3;
n = 251;
x = jacpts(n,a, b);
pass(11) = all(size(x) == [n, 1]);
[x, w, v] = jacpts(n, a, b);
pass(12) = all(size(x) == [n, 1]) && all(size(w) == [1, n]) && all(size(v) == [n, 1]);
pass(13) = abs(w*x - 5.693053859952719) < tol && abs(w*x.^2 - 5.218632704956658) < tol;
pass(14) = abs(x(37) + 0.893103435898983) < tol;
pass(15) = abs(w(37) - 1.962534523788093e-04) < tol;
pass(16) = abs(v(37) - 0.042040484981163) < tol;

% Test on [0, 10]:
[x, w, v] = jacpts(n, a, b, [0, 10]);
pass(17) = abs(w*x - 859.7954446702054) < 100*tol && abs(w*x.^2 - 7881.458242810217) < 1000*tol;
pass(18) = abs(x(38) - 0.562894092705871) < tol;
pass(19) = abs(w(38) - 0.002830847963541) < tol;
pass(20) = abs(v(38) + 0.045147639085330) < tol;

end
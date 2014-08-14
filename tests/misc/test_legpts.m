function pass = test_legpts(pref)

% Choose a tolerance:
tol = 1e-14;

% Test a small n (using REC)
n = 42;
[x] = legpts(n);
pass(1) = all(size(x) == [n, 1]);
[x, w, v] = legpts(n);
pass(2) = all(size(x) == [n, 1]) && all(size(w) == [1, n]) && all(size(v) == [n, 1]);
pass(3) = abs(w*x) < tol && abs(w*x.^2 - 2/3) < tol;
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

% Test n = 1: 
[x, w]= legpts( 1 ); 
pass(21) = (x == 0); 
pass(22) = (w == 2); 
[x, w] = legpts(1, [-10, 3]); 
pass(23) = abs( (x + 10) + (x-3) ) < tol ; % midpoint 
pass(24) = abs( w*x + 45.5 ) < tol;        % integrates x correctly on [-10,3].

% Test n = 2: 
[x, w]= legpts( 2 ); 
pass(25) = all(abs(x - [-1;1]/sqrt(3))<tol); 
pass(26) = all(abs(w - [1 1]) < tol); 
[x, w] = legpts(2, [-10, 3]); 
pass(27) = abs(sum(w) - diff([-10 3])) < tol; % integrate 1. 
pass(28) = abs( w*x + 45.5 ) < tol;            % integrates x on [-10,3].
pass(29) = abs( w*x.^2 - (342 + 1/3) ) < 20*tol;  % integrate x^2 on [-10 3]. 
pass(30) = abs( w*x.^3 + 2479.75 ) < 1000*tol;  % integrate x^3 on [-10 3]. 
end

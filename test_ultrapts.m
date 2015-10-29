function pass = test_ultrapts(pref)

% Choose a tolerance:
tol = 1e-14;

% Test a small n (using REC)
n = 42;
lambda = .3;
x = ultrapts(n, lambda);
pass(1) = all(size(x) == [n, 1]);
[x, w, v] = ultrapts(n, lambda);
pass(2) = all(size(x) == [n, 1]) && all(size(w) == [1, n]) && all(size(v) == [n, 1]);
pass(3) = abs(w*x) < tol && abs(w*x.^2 - 0.8843414686338345) < tol;
pass(4) = abs(x(37) - 9.131896381993957E-01) < tol;
pass(5) = abs(w(37) - 4.332670514309510E-02) < tol;
pass(6) = abs(v(37) - 3.115587460502451E-01) < tol;

% Test on [0, 10]:
[x, w, v] = ultrapts(n, lambda, [0, 10]);
pass(7) = abs(w*x - 30.1957169274024) < 31*tol && abs(w*x.^2 - 209.0472710358626) < 300*tol;
pass(8) = abs(x(38) - 9.704505777068543E+00) < tol;
pass(9) = abs(w(38) - 1.018229378664342E-01) < tol;
pass(10) = abs(v(38) + 2.449177929215358E-01) < tol;

% Test a larger n (using ASY without guaranteed convergence)
lambda = .8;
n = 251;
x = ultrapts(n, lambda);
pass(11) = all(size(x) == [n, 1]);
[x, w, v] = ultrapts(n, lambda);
pass(12) = all(size(x) == [n, 1]) && all(size(w) == [1, n]) && all(size(v) == [n, 1]);
pass(13) = abs(w*x) < tol && abs(w*x.^2 - 0.4744211549960596) < tol;
pass(14) = abs(x(37) + 8.958806879214126E-01) < tol;
pass(15) = abs(w(37) - 3.406945649865882E-03) < tol;
pass(16) = abs(v(37) - 2.321704534650446E-01) < tol;

% Test on [0, 10]:
[x, w, v] = ultrapts(n, lambda, [0, 10]);
pass(17) = abs(w*x - 112.1472319135050) < 100*tol && abs(w*x.^2 - 716.4962038918374) < 100*tol;
pass(18) = abs(x(38) - 5.486606034997460E-01 ) < tol;
pass(19) = abs(w(38) - 4.655102134393607E-02) < tol;
pass(20) = abs(v(38) + 2.427562206703888E-01 ) < tol;

% Test a larger n (using ASY with guaranteed convergence)
lambda = .8;
n = 251;
x = ultrapts(n, lambda,'asy',1);
pass(21) = all(size(x) == [n, 1]);
[x, w, v] = ultrapts(n, lambda,'asy',1);
pass(22) = all(size(x) == [n, 1]) && all(size(w) == [1, n]) && all(size(v) == [n, 1]);
pass(23) = abs(w*x) < tol && abs(w*x.^2 - 0.4744211549960596) < tol;
pass(24) = abs(x(37) + 8.958806879214126E-01) < tol;
pass(25) = abs(w(37) - 3.406945649865882E-03) < tol;
pass(26) = abs(v(37) - 2.321704534650446E-01) < tol;

% Test on [0, 10]:
[x, w, v] = ultrapts(n, lambda, [0, 10], 'asy', 1);
pass(27) = abs(w*x - 112.1472319135050) < 100*tol && abs(w*x.^2 - 716.4962038918374) < 100*tol;
pass(28) = abs(x(38) - 5.486606034997460E-01) < tol;
pass(29) = abs(w(38) - 4.655102134393607E-02) < tol;
pass(30) = abs(v(38) + 2.427562206703888E-01) < tol;

% Test n = 1: 
lambda = .6;
[x, w] = ultrapts(1, lambda); 
pass(31) = (x == 0); 
pass(32) = abs(w - 1.887181162535959) < tol; 
[x, w] = ultrapts(1, lambda, [-10, 3]); 
pass(33) = abs( (x + 10) + (x-3) ) < tol ; % midpoint 
pass(34) = abs( w*x + 62.42774750787926 ) < 2*tol;  

% Test n = 2: 
[x, w] = ultrapts(2, lambda); 
pass(35) = all(abs(x - [-sqrt(5)/4; sqrt(5)/4]) < tol); 
pass(36) = all(abs(w - [.5*gamma(lambda+.5)*sqrt(pi)/gamma(lambda+1),...
    .5*gamma(lambda+.5)*sqrt(pi)/gamma(lambda+1)]) < tol); 
[x, w] = ultrapts(2, lambda, [-10, 3]); 
pass(37) = abs(sum(w) - 17.83649928796550) < tol; 
pass(38) = abs( w*x + 62.42774750787926 ) < 10*tol;            
pass(39) = abs( w*x.^2 - 453.9946459389969 ) < 100*tol;  
pass(40) = abs( w*x.^3 + 3237.463968416426 ) < 100*tol;   
end

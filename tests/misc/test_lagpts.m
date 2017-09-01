


function pass = test_lagpts(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% TODO: These values were computed in V4. Should perhaps check more carefully?

% Choose a tolerance:
tol = 10*pref.chebfuneps;

% Test a small n (using GW with recur)
n = 42;
[x] = lagpts(n);
pass(1) = all(size(x) == [n, 1]);
[x, w, v] = lagpts(n);
pass(2) = all(size(x) == [n, 1]) && all(size(w) == [1, n]) && ...
    all(size(v) == [n, 1]);

pass(3) = (abs(w*x - 1) <= tol) && (abs(w*x.^2 - 2) <= tol);
pass(4) = abs(x(37) - 98.388267163326702) < 100*tol;
pass(5) = abs(w(7) - 0.055372813167092) < tol;
pass(6) = abs(v(17) - 0.002937421407003) < tol;

% Test a larger n (using RH/fast)
n = 251;
[x] = lagpts(n);
pass(7) = all(size(x) == [n, 1]);
[x, w, v] = lagpts(n);
pass(8) = all(size(x) == [n, 1]) && all(size(w) == [1, n]) && ...
    all(size(v) == [n, 1]);
pass(9) = abs(w*x - 1) < 200*tol && abs(w*x.^2 - 2) < 400*tol;
pass(10) = abs(x(37) - 13.309000189442097) < 10*tol;
pass(11) = abs(w(3) - 0.050091759039996) < 200*tol;
pass(12) = abs(v(3) - 0.214530194346947) < 200*tol;

% Test a different interval (using GW with recur)
n = 42;
[x, w, v] = lagpts(n, [1, inf]);
e = exp(1);
pass(13) = abs(w*x - 2/e) < tol && abs(w*x.^2 - 5/e) < tol;

% Put the inf on the left:
[x, w, v] = lagpts(n, [-inf, -1]);
e = exp(1);
pass(14) = abs(w*x + 2/e) < tol && abs(w*x.^2 - 5/e) < tol;

% Tests with other alpha at different n
alpha = 0.7;
n = 31; % GW/rec
[x,w] = lagpts(n,[0,inf], 'default', alpha);
exa= gamma(alpha+1) +2.5*gamma(alpha + round(n/2)+1) -31.8*gamma(alpha+n+16);
pass(15) = abs((exa - (w*x.^0 + 2.5*w*x.^round(n/2) -31.8*w*x.^(n+15)))/exa) < 10*tol;
n = 289; % RH/fast
[x,w] = lagpts(n,[0,inf], 'default', alpha);
exa= gamma(alpha+1) +2.5*gamma(alpha + 14) -31.8*gamma(alpha+56);
pass(16) = abs((exa - (w*x.^0 + 2.5*w*x.^13 -31.8*w*x.^55))/exa) < 200*tol;
n = 1681; % expl
[x,w] = lagpts(n,[0,inf], 'default', alpha);
exa= gamma(alpha+1) +2.5*gamma(alpha + 14) -31.8*gamma(alpha+29);
pass(17) = abs((exa - (w*x.^0 + 2.5*w*x.^13 -31.8*w*x.^28))/exa) < 20*tol;

n = 5001;
[x, w] = lagpts(n);
pass(18) = all(size(x) == [n,1]);
pass(19) = (abs(w*x -1) <= 100*tol);
pass(20) = (min(x) > 0.0) && ( max(x) < 4*n + 2*alpha + 2 );
pass(21) = (min(diff(x)) > 0.0) && (min(w) >= 0.0);

n = 400;
alpha = 3; % RH changes to expl
[x, w] = lagpts(n, alpha);
% High tolerance for moderate alpha
pass(22) = all(size(x) == [n,1]) && (abs(w*x -gamma(alpha+2) ) <= 1e-9);
pass(23) = (min(x) > 0.0) && ( max(x) < 4*n + 2*alpha + 2 );
pass(24) = (min(diff(x)) > 0.0) && (min(w) >= 0.0);

end


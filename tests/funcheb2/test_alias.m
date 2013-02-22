function pass = test_alias(varargin)

% [TODO]: Test aliasing a large tail (multiple aliasing).

% Set a tolerance (pref.eps doesn't matter)
tol = 100*eps;

%%
% Testing a vector of coefficients.
c0 = (1:10)';

% Padding:
c1 = funcheb.alias(c0, 11);
pass(1) = norm([0 ; c0] - c1, inf) == 0;

% Aliasing:
c2 = funcheb.alias(c0, 9);
pass(2) = norm([2, 4, 4:10]' - c2, inf) == 0;
c3 = funcheb.alias(c0, 3);
pass(3) = norm([12 25 18]' - c3, inf) == 0;

% Compare against result of evaluating on a smaller grid:
pass(4) = norm(funcheb2.chebpoly(funcheb.clenshaw(funcheb2.chebpts(9), c0)) - c2, inf) < tol;
pass(5) = norm(funcheb2.chebpoly(funcheb.clenshaw(funcheb2.chebpts(3), c0)) - c3, inf) < tol;


%%
% Testing a matrix of coefficients.
cc = [c0 c0(end:-1:1)];

% Padding:
c1 = funcheb.alias(cc, 11);
pass(6) = norm([0 0 ; cc]  - c1, inf) == 0;

% Aliasing:
c2 = funcheb.alias(cc, 9);
pass(7) = norm([2 4 4:10 ; 9 18 7:-1:1]'  - c2, inf) == 0;
c3 = funcheb.alias(cc, 3);
pass(8) = norm([12 25 18 ; 10 30 15]'  - c3, inf) == 0;

% Compare against result of evaluating on a smaller grid:
pass(9) = norm(funcheb2.chebpoly(funcheb.clenshaw(funcheb2.chebpts(9), cc)) - c2, inf) < tol;
pass(10) = norm(funcheb2.chebpoly(funcheb.clenshaw(funcheb2.chebpts(3), cc)) - c3, inf) < tol;

end

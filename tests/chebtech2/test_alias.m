% Test file for chebtech2/alias.m

function pass = test_alias(varargin)

% Set a tolerance (pref.eps doesn't matter)
tol = 100*eps;

%%
% Testing a vector of coefficients.
c0 = (10:-1:1)';

% Padding:
c1 = chebtech2.alias(c0, 11);
pass(1) = norm([c0 ; 0] - c1, inf) == 0;

% Aliasing:
c2 = chebtech2.alias(c0, 9);
pass(2) = norm([10:-1:4, 4, 2].' - c2, inf) == 0;
c3 = chebtech2.alias(c0, 3);
pass(3) = norm([18 25 12].' - c3, inf) == 0;

% Compare against result of evaluating on a smaller grid:
pass(4) = norm(chebtech2.vals2coeffs( ...
    chebtech.clenshaw(chebtech2.chebpts(9), c0)) - c2, inf) < tol;
pass(5) = norm(chebtech2.vals2coeffs( ...
    chebtech.clenshaw(chebtech2.chebpts(3), c0)) - c3, inf) < tol;


%%
% Testing a matrix of coefficients.
cc = [c0 c0(end:-1:1)];

% Padding:
c1 = chebtech2.alias(cc, 11);
pass(6) = norm([cc; 0 0]  - c1, inf) == 0;

% Aliasing:
c2 = chebtech2.alias(cc, 9);
pass(7) = norm([10:-1:4 4 2 ; 1:7 18 9]'  - c2, inf) == 0;
c3 = chebtech2.alias(cc, 3);
pass(8) = norm([18 25 12 ; 15 30 10]'  - c3, inf) == 0;

% Compare against result of evaluating on a smaller grid:
pass(9) = norm(chebtech2.vals2coeffs( ...
    chebtech.clenshaw(chebtech2.chebpts(9), cc)) - c2, inf) < tol;
pass(10) = norm(chebtech2.vals2coeffs( ...
    chebtech.clenshaw(chebtech2.chebpts(3), cc)) - c3, inf) < tol;

%%
% Test aliasing a large tail.
c0 = 1./(1:1000).^5.';
n = 17;
c1 = chebtech2.alias(c0, n);
% This should give the same result as evaluating via bary.
v0 = chebtech2.coeffs2vals(c0);
v2 = chebtech2.bary(chebtech2.chebpts(n), v0);
c2 = chebtech2.vals2coeffs(v2);
% Check in the infinity norm:
pass(11) = norm(c1 - c2, inf) < n*tol;



end

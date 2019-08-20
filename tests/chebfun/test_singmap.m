% function pass = test_singmap(pref)

% if ( nargin == 0 )
%     pref = chebfunpref();
% end

tol = 1e-13;

%% Test 1

f = @(x) exp(1+sqrt(1+x));
fcheb = chebfun(f, [-1 1], 'singmap', [.5, 0]);
F = @(x) 2*exp(1 + sqrt(1 + x)).*(-1 + sqrt(1 + x));
F = @(x) F(x) - F(-1);

x = linspace(-.9, .9);
% Test length:
pass(1) = length(f) < 100;
% Test accuracy:
pass(2) = norm(f(x) - fcheb(x), inf) < tol;
% Test sum:
pass(3) = abs(sum(fcheb) - F(1)) < tol;
% Test cumsum:
Fcheb = cumsum(fcheb);
pass(4) = norm(F(x) - Fcheb(x), inf) < tol;

% Test join:
ff = join(fcheb,fcheb);
pass(5) = norm(ff(x) - f(x), inf) < tol && norm(ff(x+2) - f(x), inf) < tol;

% Test roots:
% opts = optimset('tolx', eps, 'tolfun', eps); 
% r = fzero(@(x) f(x) - 5, 0, opts);
r = -0.628585430887966;
pass(6) = abs(roots(fcheb-5) - r) < tol;

%% Test 2

f = @(x) exp(x)+cos(x)./sqrt(x) + sqrt(1-x);
fcheb = chebfun(f, [0, 1], 'singmap', [.5, .5], 'exps', [-1 0]);
x = linspace(.1, .9);
% Test length:
pass(7) = length(f) < 100;
% Test accuracy:
pass(8) = norm(f(x) - fcheb(x), inf) < 100*tol;
% Test sum:
F1 = 4.1939969709262560649765; % Wolfram Alpha
pass(9) = abs(sum(fcheb) - F1) < 100*tol;

% Test join:
ff = join(fcheb,fcheb);
pass(10) = norm(ff(x) - f(x), inf) < 100*tol && norm(ff(x+1) - f(x), inf) < 100*tol;

% Test roots:
% opts = optimset('tolx', eps, 'tolfun', eps); 
% r = fzero(@(x) f(x) - 5, .2, opts);
r = 0.114360282839375;
pass(11) = abs(roots(fcheb-5) - r) < tol;

% end

pass
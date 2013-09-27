% A simple example BVP.

close all, clc, clear all, clear classes

% The CHEBOP:
N = chebop(@(x, u) diff(u, 2) - (u.^2), [-1, 1]);
N.lbc = @(u, v) u - 1;
% N.bc = @(x, u) sum(u);
N.bc = @(x, u) feval(u, 1) - 1;
x = chebfun('x');

rhs = x.^2;
u = N\rhs;

u = u{1}

%%
% Display:
figure(1)
plot(u, 'b'); shg
title('solutions')

%%
disp('norm of residuals:')
res = N.op(x, u) - rhs;
norm(res, 2)

disp('bc residuals:')
feval(N.lbc(u), -1)
N.bc(x, u)

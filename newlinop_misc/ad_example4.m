% A simple example BVP.

close all, clc, clear all, clear classes

pic

% The CHEBOP:
N = chebop(@(x, u) diff(u, 2) + u, [-1, 1]);
N.lbc = @(u, v) u - 1;
N.bc = @(x, u) sum(u);
x = chebfun('x');

rhs = chebmatrix({x});
u = N\rhs

poc

u = u{1}

%%
% Display:
figure(1)
plot(u, 'b'); shg
title('solutions')

%%
disp('norm of residuals:')
res = N.op(x, u) - x;
norm(res, 2)

disp('bc residuals:')
feval(N.lbc(u), -1)
N.bc(x, u)

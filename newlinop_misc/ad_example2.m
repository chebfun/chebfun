% A simple example BVP.

close all, clc, clear all, clear classes

% The CHEBOP:
N = chebop(@(x, u, v) [diff(u, 1) + v.^2 ; diff(v, 1) + sign(x).*u], [-1, 1]);
N.lbc = @(u, v) u - 1;
N.bc = @(x, u, v) sum(v);
x = chebfun('x');

uv = N\[x ; x];

%%
% Display:
u = uv{1};
v = uv{2};
figure(1)
plot(u, 'b'); hold on
plot(v, 'r'), shg, hold off
title('solutions')

%%

figure(2)
plot(diff(u, 2), 'b'); hold on
plot(diff(v), 'r'), shg, hold off
title('derivatives of solutions')

%%
disp('norm of residuals:')
res = N.op(x, u, v) -[x ; x];
norm(res{1}, 2)
norm(res{2}, 2)

disp('bc residuals:')
feval(N.lbc(u, v), -1)
N.bc(x, u, v)

% A simple example BVP.
clc, clear all, clear classes

% The 'CHEBOP':
N.op = @(x, u, v) [diff(u, 1) + v ; diff(v, 1) + sign(x).*u];
N.lbc = @(u, v) u - 1;
% N.bc = @(x, u, v) feval(v, .5);
N.bc = @(x, u, v) sum(v);

%%
% Initialise an ADCHEBFUN:
dom = [-1 0 1];
x = chebfun(@(x) x, dom);
u = adchebfun(@(u) 0*u, dom);
v = adchebfun(@(v) 0*v, dom);

u = seed(u, 1, 2);
v = seed(v, 2, 2);

%%
% Evaluate the operators to get a linearisation:
w = N.op(x, u, v);
A = vertcat(w.jacobian);
Ares = vertcat(w.func);

w = feval(N.lbc(u, v), -1);
B1 = w.jacobian;
B1res = w.func;

w = N.bc(x, u, v);
B2 = w.jacobian;
B2res = w.func;

%%
% Add BCs
A = bc(A, B1, -B1res);
A = bc(A, B2, -B2res);

% Make a CHEBMATRIX RHS:
rhs = [x ; x] - Ares;

%%
% Solve:
uv = A\rhs;

%%
% Display:
u = uv.blocks{1};
v = uv.blocks{2};
plot(u, 'b'); hold on
plot(v, 'r'), shg, hold off
title('solutions')

u(-1) - 1
% v(0.5)
sum(v)

figure
plot(diff(u, 2), 'b'); hold on
plot(diff(v), 'r'), shg, hold off
title('derivatives of solutions')


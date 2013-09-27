% A simple example BVP.

close all, clc, clear all, clear classes

% The CHEBOP:
N = chebop(@(x, u, v) [diff(u, 1) + v.^2 ; diff(v, 1) + sign(x).*u], [-1, 1]);
% N.lbc = @(u, v) u - 1;
% N.bc = @(x, u, v) sum(v);
% N.lbc = @(u, v) [u - 1 ; v - 1];
N.bc = @(x, u, v) [feval(u,-1)-1 ; sum(v)];


N = chebop(@(x, u, v) [diff(u, 2) + v.^2 ; diff(v, 2) + sign(x).*u], [-1, 1]);
% N.lbc = @(u, v) u - 1;
% N.bc = @(x, u, v) sum(v);
N.rbc = @(u, v) [u - 1 ; v - 1];
N.bc = @(x, u, v) [feval(u,-1)-1 ; sum(v)];

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
plot(diff(u, 3), 'b'); hold on
plot(diff(v, 2), 'r'), shg, hold off
title('derivatives of solutions')

%%
disp('norm of residuals:')
res = N.op(x, u, v) -[x ; x];
norm(res{1}, 2)
norm(res{2}, 2)

disp('bc residuals:')
if ( ~isempty(N.lbc) )
    feval(N.lbc(u, v), N.domain(1))
end
if ( ~isempty(N.rbc) )
    feval(N.rbc(u, v), N.domain(end))
end
if ( ~isempty(N.bc) )
    N.bc(x, u, v)
end

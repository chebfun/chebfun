% A simple example BVP.

close all, clc, clear all, clear classes

dom = [0 1];

% The CHEBOP:
N = chebop(@(x, u) diff(u, 2) + sin(u), dom);
N.lbc = @(u, v) u - 2;
% N.bc = @(x, u) sum(u);
N.bc = @(x, u) feval(u, dom(2)) - 2;
x = chebfun('x', dom);
% N.init = 2+0*x;
N.init = 2*cos(2*pi*x/10);

rhs = 0*x;
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
norm(res, inf)

disp('bc residuals:')
if ( ~isempty(N.lbc) )
    feval(N.lbc(u), N.domain(1))
end
if ( ~isempty(N.rbc) )
    feval(N.rbc(u), N.domain(end))
end
if ( ~isempty(N.bc) )
    N.bc(x, u)
end

%% First order demo
% AB, 2014/07/24
%

%% Colour coded variables when plotting expression trees:
u = treeVar([1 0 0]);
v = treeVar([0 1 0]);
w = treeVar([0 0 1]);
f = diff(v) + v - 28*u + u.*w; 
treeVar.plotTree(f.tree)

%% First order initial value problem
dom = [0, 1];
N = chebop(@(u) diff(u) + u, dom);
N.lbc = 1;
u = N\0
plot(u), shg
fprintf('Norm of residual: %4.2e.\n', norm(N(u)))
fprintf('Residual of initial condition: %4.2e.\n\n', ...
    norm(feval(N.lbc(u), dom(1))))

%% First order initial value problem, with different coefficient on u'
dom = [0, 1];
N = chebop(@(u) .1*diff(u) + u, dom);
N.lbc = 1;
u = N\0
plot(u), shg
fprintf('Norm of residual: %4.2e.\n', norm(N(u)))
fprintf('Residual of initial condition: %4.2e.\n\n', ...
    norm(feval(N.lbc(u), dom(1))))

%% Problem with a variable coefficient
dom = [0, 2];
N = chebop(@(x, u) .1*diff(u) + sin(4*x).*u, dom);
N.lbc = 1;
u = N\0
plot(u), shg
fprintf('Norm of residual: %4.2e.\n', norm(N(u)))
fprintf('Residual of initial condition: %4.2e.\n\n', ...
    norm(feval(N.lbc(u), dom(1))))

%% Problem with a variable coefficient and an affine part
dom = [0, 2];
N = chebop(@(x, u) .1*diff(u) + sin(4*x).*u + .05*cos(x), dom);
N.lbc = 1;
u = N\0
plot(u), shg
fprintf('Norm of residual: %4.2e.\n', norm(N(u)))
fprintf('Residual of initial condition: %4.2e.\n\n', ...
    norm(feval(N.lbc(u), dom(1))))

%% Final value problem:
dom = [0, 1];
N = chebop(@(u) 0.5*diff(u) + u, dom);
N.rbc = 1;
u = N\0
plot(u), shg
fprintf('Norm of residual: %4.2e.\n', norm(N(u)))
fprintf('Residual of final condition: %4.2e.\n\n', ...
    norm(feval(N.rbc(u), dom(end))))

%% Second order initial value problem
dom = [0, 1];
N = chebop(@(x,u) diff(u, 2) + u, dom);
N.lbc = @(u) [u-1; diff(u)-1];
u = N\0
plot(u), shg
fprintf('Norm of residual: %4.2e.\n', norm(N(u)))
fprintf('Residual of initial condition: %4.2e.\n\n', ...
    norm(feval(N.lbc(u), dom(1))))

%% Problem with a variable coefficient
dom = [0, 2];
N = chebop(@(x, u) .1*diff(u, 2) + sin(4*x).*u, dom);
N.lbc = @(u) [u-1; diff(u)];
u = N\0
plot(u), shg
fprintf('Norm of residual: %4.2e.\n', norm(N(u)))
fprintf('Residual of initial condition: %4.2e.\n\n', ...
    norm(feval(N.lbc(u), dom(1))))


%% Problem with a variable coefficient and an affine part
dom = [0, 2];
N = chebop(@(x, u) .1*diff(u, 2) + sin(4*x).*u + .05*cos(x), dom);
N.lbc = @(u) [u-1; diff(u)];
u = N\0
plot(u), shg
fprintf('Norm of residual: %4.2e.\n', norm(N(u)))
fprintf('Residual of initial condition: %4.2e.\n\n', ...
    norm(feval(N.lbc(u), dom(1))))

%% van der Pol
mu = 10;
dom = [0, 50];
vdpFun = @(y) diff(y, 2) - mu.*(1-y.^2).*diff(y) + y;
N = chebop(vdpFun, dom);
N.lbc = @(u) [u-2; diff(u)];
y = N\0;
plot(y), shg
% Show where breakpoints got introduced
hold on, plot(y.domain, y(y.domain),'k*')

%% Second order final value problem
dom = [0, 1];
N = chebop(@(u) 0.5*diff(u, 2) + u, dom);
N.rbc = @(u) [u-2; diff(u)];
u = N\0
plot(u)
fprintf('Norm of residual: %4.2e.\n', norm(N(u)))
fprintf('Residual of final condition: %4.2e.\n\n', ...
    norm(feval(N.rbc(u), dom(end))))
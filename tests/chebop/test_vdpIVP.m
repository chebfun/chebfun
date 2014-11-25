function pass = test_vdpIVP(~)

%% Setup
vdpFun = @(mu) @(u) diff(u, 2) - mu.*(1-u.^2).*diff(u) + u;
dom = [0 10];
bcFun =  @(u) [u-2; diff(u)];
opts = odeset('abstol', 1e5*eps, 'reltol', 100*eps);
tol = 1e-10;

%% Solve with CHEBOPs, compare with standard MATLAB, mu = 1
N = chebop(vdpFun(1), dom);
N.lbc = bcFun;
u = N\0;
uend = u(dom(end));
sol =ode45(@vdp1, dom, [2 0], opts);  
yend = deval(sol, dom(end)); yend = yend(1);
err(1) = abs(uend - yend);

%% Solve with CHEBOPs, compare with standard MATLAB, mu = 1000
N = chebop(vdpFun(1000), dom);
N.lbc = bcFun;
u = N\0;
uend = u(dom(end));
sol =ode45(@vdp1000, dom, [2 0], opts);  
yend = deval(sol, dom(end)); yend = yend(1);
err(2) = abs(uend - yend);

%% Happy?
pass = err < tol;

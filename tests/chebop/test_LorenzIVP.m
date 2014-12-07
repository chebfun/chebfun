function pass = test_LorenzIVP(~)
%TEST_LORENZIVP   Solve the Lorenz IVP using both CHEBOPs and standard MATLAB

%% Setup
dom = [0 5];
tol = 1e-14;

%% With CHEBOP
N = chebop(@(t,u,v,w) [diff(u) - 10*(v - u);
    diff(v) - u.*(28 - w) + v;
    diff(w) - u.*v + (8/3)*w], dom);
N.lbc = @(u,v,w) [w - 20 ; v + 15; u + 14];
uvw = N\[0;0;0];
u = uvw{1}; v = uvw{2}; w = uvw{3};

%% Using standard MATLAB
odeFun = @(t,y) [10*(y(2) - y(1));
    y(1).*(28-y(3)) - y(2);
    y(1).*y(2) - (8/3)*y(3)];
opts = odeset('abstol', 1e5*eps, 'reltol', 100*eps);
sol = ode113(odeFun, dom, [-14, -15, 20], opts);
solEnd = deval(sol, dom(end));

%% Compare results:
err = [u(end); v(end); w(end)] - solEnd;

%% Happy?
pass = norm(err) < tol;
end
% Test file for SPIN3:

function pass = test_spin3()

tol = 1e-2;

%% GL:

% Solve with DT and DT/2:
S = spinop3('GL'); S.tspan = S.tspan/10;
N = 32; dt = 1e-1;
s = rng; u = spin3(S, N, dt, 'plot', 'off');
rng(s); v = spin3(S, N, dt/2, 'plot', 'off');

% Compare:
dom = S.domain;
pts = trigtech.tensorGrid([20 20 20], dom);
xx = pts{1}; yy = pts{2}; zz = pts{3};
scale = max(max(max(abs(v(xx,yy,zz)))));
pass(1) = max(max(max(abs(u(xx,yy,zz) - v(xx,yy,zz)))))/scale < tol;

end
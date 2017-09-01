% Test file for SPIN2:

function pass = test_spin2()

tol = 1e-4;

%% GL:

% Solve with DT and DT/2:
S = spinop2('GL'); S.tspan = S.tspan/10;
N = 128; dt = 1e-1;
s = rng; u = spin2(S, N, dt, 'plot', 'off');
rng(s); v = spin2(S, N, dt/2, 'plot', 'off');

% Compare:
dom = S.domain;
[xx, yy] = meshgrid(linspace(dom(1), dom(2), 50));
scale = max(max(abs(v(xx,yy))));
pass(1) = max(max(abs(u(xx,yy) - v(xx,yy))))/scale < tol;

end
% Test file for SPIN:

function pass = test_spin()

tol = 1e-6;

%% AC:

% Solve with DT and DT/2:
S = spinop('AC');
N = 256; dt = 1e-1;
u = spin(S, N, dt, 'plot', 'off');
v = spin(S, N, dt/2, 'plot', 'off');

% Compare:
pass(1) = norm(u - v, inf)/norm(v, inf) < tol;

%% CH:

% Solve with DT and DT/2:
S = spinop('CH');
N = 256; dt = 2e-2;
u = spin(S, N, dt, 'plot', 'off');
v = spin(S, N, dt/2, 'plot', 'off');

% Compare:
pass(2) = norm(u - v, inf)/norm(v, inf) < tol;

end
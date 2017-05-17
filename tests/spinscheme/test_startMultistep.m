% Test file for SPINCHEME/STARTMULTISTEP:

function pass = test_startMultistep()

tol = 1e-10;

%% 1D:

% Solve 1D KDV equation with a one-step method (ETDRK4) and with a multistep
% method (PECEC736) which needs to do 5 initialization time-steps:
S = spinop('kdv');
dt = 1e-7;
N = 256;
S.tspan = [0 10*dt];

% Solve with PECEC736:
umulti = spin(S, N, dt, 'plot', 'off', 'scheme', 'pecec736');

% Solve with ETDRK4: 
u = spin(S, N, dt, 'plot', 'off', 'scheme', 'etdrk4');

% Compare:
pass(1) = norm(u - umulti, inf)/norm(u, inf) < tol;

%% 2D:

% Solve 2D GS equation with a one-step method (ETDRK4) and with a multistep
% method (PECEC736) which needs to do 5 initialization time-steps:
S = spinop2('gs2');
dt = 1e-2;
N = 64;
S.tspan = [0 10*dt];

% Solve with ETDRK4 from 0 to 5*dt:
u = spin2(S, N, dt, 'plot', 'off', 'scheme', 'etdrk4');

% Solve with PECEC736 from 0 to 5*dt:
umulti = spin2(S, N, dt, 'plot', 'off', 'scheme', 'pecec736');

% Compare:
pass(2) = norm(u-umulti, inf)/norm(u, inf) < tol;

end
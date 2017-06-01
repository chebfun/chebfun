% Test file for SPINCHEME/STARTMULTISTEP:

function pass = test_startMultistep()

tol = 1e-10;

%% 1D:

% Solve 1D KDV equation with a one-step method (ETDRK4) and with a multistep
% method (PECEC736) which needs to do 5 initialization time-steps:
S = spinop('KDV');
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
S = spinop2('GS');
dt = 1e-2;
N = 64;
S.tspan = [0 10*dt];

% Solve with ETDRK4:
u = spin2(S, N, dt, 'plot', 'off', 'scheme', 'etdrk4');

% Solve with PECEC736:
umulti = spin2(S, N, dt, 'plot', 'off', 'scheme', 'pecec736');

% Compare:
pass(2) = norm(u-umulti, inf)/norm(u, inf) < tol;

%% SPHERE:

% Solve the AC equation with a one-step method (LIRK4) and with a multistep
% method (IMEXBDF4) which needs to do 3 initialization time-steps:
S = spinopsphere('AC');
dt = 5e-4;
N = 128;
S.tspan = [0 10*dt];

% Solve with LIRK4:
u = spinsphere(S, N, dt, 'plot', 'off', 'scheme', 'lirk4');

% Solve with IMEXBDF4:
umulti = spinsphere(S, N, dt, 'plot', 'off', 'scheme', 'imexbdf4');

% Compare:
pass(3) = norm(u-umulti, inf)/norm(u, inf) < tol;

end
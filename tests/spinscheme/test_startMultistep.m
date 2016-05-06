% Test file for SPINCHEME/STARTMULTISTEP:

function pass = test_startMultistep()

tol = 1e-10;

%% 1D:

% Solve 1D KDV equation with a one-step method (ETDRK4) and with a multistep
% method (PECEC736) which needs to do 5 initialization time-steps:
dt = 1e-6;
pref = spinpref('dt', dt, 'N', 256);
pref.plot = 'off';

% Solve with ETDRK4 from 0 to 5*dt:
pref.scheme = 'etdrk4';
u = spin('kdv', pref, [0 5*dt]);

% Solve with PECEC736 from 0 to 5*dt:
pref.scheme = 'pecec736';
umulti = spin('kdv', pref, [0 5*dt]);

% Compare:
pass(1) = norm(u-umulti) < tol;

%% 2D:

% Solve 2D GS equation with a one-step method (ETDRK4) and with a multistep
% method (PECEC736) which needs to do 5 initialization time-steps:
dt = 1e-2;
pref = spinpref2('dt', dt, 'N', 64);
pref.plot = 'off';

% Solve with ETDRK4 from 0 to 5*dt:
pref.scheme = 'etdrk4';
u = spin2('gs2', pref, [0 5*dt]);

% Solve with PECEC736 from 0 to 5*dt:
pref.scheme = 'pecec736';
umulti = spin2('gs2', pref, [0 5*dt]);

% Compare:
pass(2) = norm(u-umulti) < tol;

end

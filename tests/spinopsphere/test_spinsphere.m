% Test file for SPINSPHERE:

function pass = test_spinsphere()

tol = 1e-4;

%% GL:

% Create a SPINOPSPHERE:
S = spinopsphere('gl');
dt = 1e-2;
S.tspan = [0 100*dt];
S.init = spherefun(@(x,y,z) cos(x) + cos(y) + cos(z));

% Solve with dt:
u = spinsphere(S, 32, dt, 'plot', 'off');

% Solve with dt/2 and check error is small:
uexact = spinsphere(S, 32, dt/2, 'plot', 'off');

% Compare:
tic, norm(u - uexact), toc
pass(1) = norm(u - uexact) < tol;

end
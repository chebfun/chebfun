
% Test file for SPIN:

function pass = test_spin()

tol = 1e-6;

%% KDV:

% Initial condition:
dom = [-pi pi];
A = 25; 
B = 16;
u0 = @(x) 3*A^2*sech(.5*A*x).^2 + 3*B^2*sech(.5*B*(x-1)).^2;
u0 = chebfun(u0, dom);

% Solve up to T=4e-3:
T = 4e-3;
pref = spinpref('dt', 1e-6, 'N', 256, 'plot', 'off');
u = spin('kdv', [0 T], u0, pref);

% Exact solution at T=4e-3:
phiA = 1/A*log((A + B)^2/(A - B)^2);
phiB = 1/B*log((A + B)^2/(A - B)^2);
uexact = chebfun(@(x) 3*A^2*sech(.5*A*(x - A^2*T - phiA)).^2 + ...
    3*B^2*sech(.5*B*(x - 1 - B^2*T + phiB)).^2, dom);

% Compare:
pass(1) = (norm(u - uexact, inf)/norm(uexact, inf)) < 100*tol;

%% NLS:

% Initial condition:
dom = [-pi pi];
b = 1; 
a = 2;
u0 = @(x) (2*b^2./(2 - sqrt(2)*sqrt(2-b^2)*cos(a*b*x)) - 1)*a;
u0 = chebfun(u0, dom, 'trig');

% Solve up to T=1:
T = 1;
pref = spinpref('dt', 1e-3, 'N', 256, 'plot', 'off');
u = spin('NLS', [0 T], u0, pref);
 
% Exact solution at T=1:
theta = a^2*b*sqrt(2 - b^2)*T;
uexact = @(x) ((2*b^2*cosh(theta) + 2i*b*sqrt(2-b^2)*sinh(theta))./...
    (2*cosh(theta) - sqrt(2)*sqrt(2-b^2)*cos(a*b*x)) - 1).*(a*exp(1i*a^2*T));
uexact = chebfun(uexact, dom, 'trig');

% Compare:
pass(2) = (norm(u - uexact, inf)/norm(uexact, inf)) < tol;

%% AC:

% Solve AC with dt up to T=10:
T = 10;
dt = 2e-2;
pref = spinpref('dt', dt, 'N', 256, 'plot', 'off');
u = spin('ac', pref, [0 T]);

% Solve with dt/2 and check error is small:
pref.dt = dt/2;
uexact = spin('ac', pref, [0 T]);

% Compare:
pass(3) = (norm(u - uexact, inf)/norm(uexact, inf)) < tol;

%% AC (BIS): TEST PARSING

% Same equation and same preferences, but with different syntax:
v = spin('ac', [0 T], 'dt', dt, 'N', 256, 'plot', 'off');

% Compare: (u and v should be the same);
pass(4) = norm(u - v, inf) == 0;

end

% Test file for SPIN:

function pass = test_spin()

tol = 1e-6;

%% KDV:

% Create a SPINOP:
S = spinop('kdv');
T = 4e-3;
S.tspan = [0 T];
A = 25; 
B = 16;
dom = S.domain;
S.init = chebfun(@(x) 3*A^2*sech(.5*A*x).^2 + 3*B^2*sech(.5*B*(x-1)).^2, dom);

% Solve:
u = spin(S, 256, 1e-6, 'plot', 'off');

% Exact solution:
phiA = 1/A*log((A + B)^2/(A - B)^2);
phiB = 1/B*log((A + B)^2/(A - B)^2);
uexact = chebfun(@(x) 3*A^2*sech(.5*A*(x - A^2*T - phiA)).^2 + ...
    3*B^2*sech(.5*B*(x - 1 - B^2*T + phiB)).^2, dom);

% Compare:
pass(1) = (norm(u - uexact, inf)/norm(uexact, inf)) < 100*tol;

%% NLS:

% Create a SPINOP:
S = spinop('nls');
T = 1;
S.tspan = [0 1];

% Solve:
u = spin(S, 256, 1e-3, 'plot', 'off');
 
% Exact solution:
b = 1; 
a = 2;
theta = a^2*b*sqrt(2 - b^2)*T;
uexact = @(x) ((2*b^2*cosh(theta) + 2i*b*sqrt(2-b^2)*sinh(theta))./...
    (2*cosh(theta) - sqrt(2)*sqrt(2-b^2)*cos(a*b*x)) - 1).*(a*exp(1i*a^2*T));
uexact = chebfun(uexact, S.domain, 'trig');

% Compare:
pass(2) = (norm(u - uexact, inf)/norm(uexact, inf)) < tol;

%% AC:

% Create a SPINOP:
S = spinop('ac');
T = 10;
S.tspan = [0 T];

% Solve with dt:
dt = 1e-2;
u = spin(S, 256, dt, 'plot', 'off');

% Solve with dt/2 and check error is small:
uexact = spin(S, 256, dt/2, 'plot', 'off');

% Compare:
pass(3) = (norm(u - uexact, inf)/norm(uexact, inf)) < tol;

end
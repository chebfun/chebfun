% Test file for SPIN2:

function pass = test_spin2()

tol = 1e-4;

%% GL2:

% Create a SPINOP2:
S = spinop2('gl2');
dom = S.domain;
T = 20;
S.tspan = [0 T];
S.init = chebfun2(@(x,y) exp(-.1*((x-50).^2 + (y-50).^2)), dom, 'trig');

% Solve with dt:
dt = 1e-1;
u = spin2(S, 64, dt, 'plot', 'off');

% Solve with dt/2 and check error is small:
uexact = spin2(S, 64, dt/2, 'plot', 'off');

% Compare:
[xx, yy] = meshgrid(linspace(0, 100, 50));
scale = max(max(abs(uexact(xx,yy))));
pass(1) = max(max(abs(u(xx,yy) - uexact(xx,yy))))/scale < tol;

%% GS2:

% Create a SPINOP2:
S = spinop2('gs2');
T = 100;
S.tspan = [0 T];

% Solve with dt:
dt = 1;
u = spin2(S, 64, dt, 'plot', 'off');

% Solve with dt/2 and check error is small:
uexact = spin2(S, 64, dt/2, 'plot', 'off');

% Compare:
[xx, yy] = meshgrid(linspace(0, 100, 50));
scale = max(max(abs(uexact{1}(xx,yy))));
pass(2) = max(max(abs(u{1}(xx,yy) - uexact{1}(xx,yy))))/scale < tol;

end
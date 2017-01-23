% Test file for SPIN3:

function pass = test_spin3()

tol = 1e-4;

%% GL3:

% Create a SPINOP3:
S = spinop3('gl3');
dom = S.domain;
T = 20;
S.tspan = [0 T];
u0 = chebfun3(@(x,y,z) exp(-((x-50).^2 + (y-50).^2 + (z-50).^2)), dom, 'trig');
S.init = u0;

% Solve with dt:
dt = 1e-1;
u = spin3(S, 20, dt, 'plot', 'off');

% Solve with dt/2 and check error is small:
uexact = spin3(S, 20, dt/2, 'plot', 'off');

% Compare:
pts = trigtech.tensorGrid([20 20 20], dom);
xx = pts{1};
yy = pts{2};
zz = pts{3};
scale = max(max(max(abs(uexact(xx,yy,zz)))));
pass(1) = max(max(max(abs(u(xx,yy,zz) - uexact(xx,yy,zz)))))/scale < tol;

%% GS3:

% Create a SPINOP3:
S = spinop3('gs3');
T = 100;
S.tspan = [0 T];

% Solve with dt:
dt = 1;
u = spin3(S, 20, dt, 'plot', 'off');

% Solve with dt/2 and check error is small:
uexact = spin3(S, 20, dt/2, 'plot', 'off');

% Compare:
pts = trigtech.tensorGrid([20 20 20], dom);
xx = pts{1};
yy = pts{2};
zz = pts{3};
scale = max(max(max(abs(uexact{1}(xx,yy,zz)))));
pass(2) = max(max(max(abs(u{1}(xx,yy,zz) - uexact{1}(xx,yy,zz)))))/scale < tol;

end
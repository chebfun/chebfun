% Test file for SPIN2:

function pass = test_spin3()

tol = 1e-4;

%% G3:

% Solve GL3 with dt up to T=20;
T = 20;
dt = 1e-1;
dom = [0 100 0 100 0 100];
u0 = chebfun3(@(x,y,z) exp(-((x-50).^2 + (y-50).^2 + (z-50).^2)), dom, 'trig');
pref = spinpref3('dt', dt, 'N', 20, 'plot', 'off');
u = spin3('gl3', pref, u0, [0 T]);

% Solve with dt/2 and check error is small:
pref.dt = dt/2;
uexact = spin3('gl3', pref, u0, [0 T]);

% Compare:
pts = trigtech.tensorGrid([20 20 20], dom);
xx = pts{1};
yy = pts{2};
zz = pts{3};
scale = max(max(max(abs(uexact(xx,yy,zz)))));
pass(1) = max(max(max(abs(u(xx,yy,zz) - uexact(xx,yy,zz)))))/scale < tol;

%% GS3:

% Solve GS3 with dt up to T=100:
T = 100;
dt = 1;
pref = spinpref3('dt', dt, 'N', 20, 'plot', 'off');
u = spin3('gs3', pref, [0 T]);

% Solve with dt/2 and check error is small:
pref.dt = dt/2;
uexact = spin3('gs3', pref, [0 T]);

% Compare:
pts = trigtech.tensorGrid([20 20 20], dom);
xx = pts{1};
yy = pts{2};
zz = pts{3};
scale = max(max(max(abs(uexact{1}(xx,yy,zz)))));
pass(2) = max(max(max(abs(u{1}(xx,yy,zz) - uexact{1}(xx,yy,zz)))))/scale < tol;

%% GS3 (BIS): TEST PARSING

% Same equation and same preferences, but with different syntax:
v = spin3('gs3', [0 T], 'dt', dt, 'N', 20, 'plot', 'off');

% Compare: (u and v should be the same);
pass(3) = norm(u - v, inf) == 0;

end

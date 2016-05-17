% Test file for SPIN2:

function pass = test_spin2()

tol = 1e-4;

%% GL2:

% Solve GL2 with dt up to T=20;
T = 20;
dt = 1e-1;
dom = [0 100 0 100];
u0 = chebfun2(@(x,y) exp(-.1*((x-50).^2 + (y-50).^2)), dom, 'trig');
pref = spinpref2('dt', dt, 'N', 64, 'plot', 'off');
u = spin2('gl2', pref, u0, [0 T]);

% Solve with dt/2 and check error is small:
pref.dt = dt/2;
uexact = spin2('gl2', pref, u0, [0 T]);

% Compare:
[xx, yy] = meshgrid(linspace(0, 100, 50));
scale = max(max(abs(uexact(xx,yy))));
pass(1) = max(max(abs(u(xx,yy) - uexact(xx,yy))))/scale < tol;

%% GS2:

% Solve GS2 with dt up to T=100:
T = 100;
dt = 1;
pref = spinpref2('dt', dt, 'N', 64, 'plot', 'off');
u = spin2('gs2', pref, [0 T]);

% Solve with dt/2 and check error is small:
pref.dt = dt/2;
uexact = spin2('gs2', pref, [0 T]);

% Compare:
[xx, yy] = meshgrid(linspace(0, 100, 50));
scale = max(max(abs(uexact{1}(xx,yy))));
pass(2) = max(max(abs(u{1}(xx,yy) - uexact{1}(xx,yy))))/scale < tol;

%% GS2 (BIS): TEST PARSING

% Same equation and same preferences, but with different syntax:
v = spin2('gs2', [0 T], 'dt', dt, 'N', 64, 'plot', 'off');

% Compare: (u and v should be the same);
pass(3) = norm(u - v, inf) == 0;

end

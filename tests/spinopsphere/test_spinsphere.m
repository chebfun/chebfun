% Test file for SPINSPHERE:

function pass = test_spinsphere()

tol = 1e-2;

%% AC:

% Solve with DT and DT/2:
S = spinopsphere('AC'); S.tspan = S.tspan/10;
N = 128; dt = 1e-1;
u = spinsphere(S, N, dt, 'plot', 'off');
v = spinsphere(S, N, dt/2, 'plot', 'off');

% Compare:
dom = S.domain;
[xx, yy] = meshgrid(linspace(dom(1), dom(2), 50));
scale = max(max(abs(v(xx,yy))));
pass(1) = max(max(abs(u(xx,yy) - v(xx,yy))))/scale < tol;

%% GL:

% Solve with DT and DT/2:
S = spinopsphere('GL'); S.tspan = S.tspan/10; 
S.init =  spherefun(@(x,y,z) cos(cosh(x.*z)-y));
N = 128; dt = 1e-1;
u = spinsphere(S, N, dt, 'plot', 'off');
v = spinsphere(S, N, dt/2, 'plot', 'off');

% Compare:
dom = S.domain;
[xx, yy] = meshgrid(linspace(dom(1), dom(2), 50));
scale = max(max(abs(v(xx,yy))));
pass(2) = max(max(abs(u(xx,yy) - v(xx,yy))))/scale < tol;

end

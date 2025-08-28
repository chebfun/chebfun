function pass = test_ivp(varargin)

% This test solves the Van der Pol ODE in CHEBFUN using ode15s and ode45. It
% checks the solution with Matlab's inbuilt ode15s and ode45 solvers.

% Rodrigo Platte Jan 2009, modified by Asgeir Birkisson, modified by Nick Hale.

%% Test using Van der Pol:

% Test ode15s Using default tolerances (RelTol = 1e-3)
[t, y] = chebfun.ode15s(@vdp1, [0, 5], [2 ; 0]); % CHEBFUN solution
[tm, ym] = ode15s(@vdp1, [0, 5], [2 ; 0]);       % Matlab's solution
pass(1) = max(max(abs(ym - feval(y, tm)))) < 2e-2;

% Test ode45 Using default tolerances (RelTol = 1e-3)
[t, y] = chebfun.ode45(@vdp1, [0, 5], [2 ; 0]);  % CHEBFUN solution
[tm, ym] = ode45(@vdp1, [0, 5], [2 ; 0]);        % Matlab's solution
pass(2) = max(max(abs(ym - feval(y, tm)))) < 1e-2;

% Test ode113 with different tolerance (RelTol = 1e-6)
opts = odeset('RelTol', 1e-6);
[t, y] = chebfun.ode113(@vdp1, [0, 20], [2 ; 0], opts); % CHEBFUN solution
[tm, ym] = ode113(@vdp1, [0, 20], [2 ; 0], opts);       % Matlab's solution
pass(3) = max(max(abs(ym - feval(y,tm)))) < 1e-5;

if ( verLessThan('matlab', '9.11') )
    pass(4:5) = true;
else
    % Test ode78 with different tolerance (RelTol = 1e-6)
    opts = odeset('RelTol', 1e-6);
    [t, y] = chebfun.ode78(@vdp1, [0, 20], [2 ; 0], opts); % CHEBFUN solution
    [tm, ym] = ode78(@vdp1, [0, 20], [2 ; 0], opts);       % Matlab's solution
    pass(4) = max(max(abs(ym - feval(y,tm)))) < 1e-5;

    % Test ode89 with different tolerance (RelTol = 1e-6)
    opts = odeset('RelTol', 1e-6);
    [t, y] = chebfun.ode89(@vdp1, [0, 20], [2 ; 0], opts); % CHEBFUN solution
    [tm, ym] = ode89(@vdp1, [0, 20], [2 ; 0], opts);       % Matlab's solution
    pass(5) = max(max(abs(ym - feval(y,tm)))) < 1e-5;
end

%% Test some trivial complex-valued IVPs:

f = @(x, u) 1i*u;
d = [0, 1];
soln = exp(1i);

% Test ode15s:
y = chebfun.ode15s(f, d, 1);
pass(6) = abs(y(1) - soln) < 2e-2;

% Test ode45:
y = chebfun.ode45(f, d, 1);
pass(7) = abs(y(1) - soln) < 2e-2;

% Test ode113:
y = chebfun.ode113(f, d, 1);
pass(8) = abs(y(1) - soln) < 2e-2;

if ( verLessThan('matlab', '9.11') )
    pass(9:10) = true;
else
    % Test ode78:
    y = chebfun.ode78(f, d, 1);
    pass(9) = abs(y(1) - soln) < 2e-2;

    % Test ode89:
    y = chebfun.ode89(f, d, 1);
    pass(10) = abs(y(1) - soln) < 2e-2;
end

end


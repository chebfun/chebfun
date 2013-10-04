function pass = test_ode(pref)

if ( nargin == 0 )
    pref = chebfun.pref();
end

tol = 1e-3;

%% IVP:

% Solve Van der Pol problem:
r_true = 3.503054040938623;

y = chebfun.ode113(@vdp1, [0, 20], [2 ; 0]); 
r1 = roots( y(:, 1) - 1 );
pass(1) = abs(r1 - r_true) < tol;


y = chebfun.ode15s(@vdp1, [0, 20], [2 ; 0]); 
r3 = roots( y(:, 1) - 1 );
pass(2) = abs(r3 - r_true) < tol;

y = chebfun.ode45(@vdp1, [0, 20], [2 ; 0]); 
r2 = roots( y(:, 1) - 1 );
pass(3) = abs(r2 - r_true) < tol;

% Test complex values:
f = @(x, u) 1i*u;
d = [0, 1];

y = chebfun.ode113(f, d, 1);
pass(4) = abs(y(1) - exp(1i)) < 2e-2;

y = chebfun.ode15s(f, d, 1);
pass(5) = abs(y(1) - exp(1i)) < 2e-2;

y = chebfun.ode45(f, d, 1);
pass(6) = abs(y(1) - exp(1i)) < 2e-2;

%% BVP:

y2_true = [1.879465902962025, -0.860039028451148];
d = [0, 4];
y0 = chebfun(@(x) [1+0*x, 0*x], d);

y = bvp4c(@twoode, @twobc, y0);
y2 = feval(y, 2);
pass(7) = norm(y2 - y2_true, inf) < tol;

y = bvp5c(@twoode, @twobc, y0);
y2 = feval(y, 2);
pass(8) = norm(y2 - y2_true, inf) < tol;

% Set new tolerance:
y0 = chebfun([0, 0], [0, 4]);
solinit = bvpinit([0, 1, 2, 3, 4], [1, 0]); 
opts = odeset('RelTol', 1e-6);
y = bvp5c(@twoode, @twobc, y0, opts);         % Chebfun solution
sol = bvp5c(@twoode, @twobc, solinit, opts);  % Matlab's solution
pass(9) = max(max(abs(sol.y' - feval(y, sol.x')))) < 1e-4; 


%% A problem with a parameter:

p_true = 17.096591689705100;
sol1_true = [-0.703689352093852, 1.033088348257337];

% Problem parameter, shared with nested functions.
q = 5;
lambda = 15;

% Derivative function. q is provided by the outer function.
mat4ode = @(x, y, lambda) [ y(2) ; -(lambda - 2*q*cos(2*x))*y(1) ];
% Boundary conditions. lambda is a required argument.
mat4bc = @(ya, yb,  lambda) [ ya(2) ; yb(2) ; ya(1)-1 ];
% Auxiliary function -- initial guess for the solution
mat4init = @(x) [ cos(4*x) , -4*sin(4*x) ];

solinit = bvpinit(linspace(0,pi,10), mat4init, lambda);
solinit = chebfun(mat4init, [0, pi]);

[sol, p] = bvp4c(mat4ode, mat4bc, solinit, lambda);
sol1 = feval(sol, 1);
pass(10) = norm(sol1 - sol1_true, inf) < tol;
pass(11) = abs(p - p_true) < tol;

[sol, p] = bvp5c(mat4ode, mat4bc, solinit, lambda);
sol1 = feval(sol, 1);
pass(10) = norm(sol1 - sol1_true, inf) < tol;
pass(11) = abs(p - p_true) < tol;

end







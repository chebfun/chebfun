function pass = test_bvp5c(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

%% Test using TWOODE():
d = [0, 4];
y0 = chebfun([1, 0], d, pref);
solinit = bvpinit([0, 1, 2, 3, 4], [1, 0]); 
% Test bvp5c using default tolerance (RelTol = 1e-3)
y = bvp5c(@twoode, @twobc, y0);         % Chebfun solution
sol = bvp5c(@twoode, @twobc, solinit);  % Matlab's solution
pass(1) = max(max(abs(sol.y' - feval(y,sol.x')))) < 2e-2;


%% Test using MAT4BC(): (A problem with a parameter)

tol = 1e-4;

p_true = 17.096591689705100;
sol1_true = [-0.703689352093852, 1.033088348257337];

% Problem parameter, shared with nested functions.
q = 5;
lambda = 15;

% Derivative function. q is provided by the outer function.
mat4ode = @(x, y, lambda) [ y(2) ; -(lambda - 2*q*cos(2*x))*y(1) ];
% Boundary conditions. lambda is a required argument.
mat4bc = @(ya, yb, lambda) [ ya(2) ; yb(2) ; ya(1)-1 ];
% Auxiliary function -- initial guess for the solution
mat4init = @(x) [ cos(4*x) , -4*sin(4*x) ];

opts = odeset('AbsTol', 1e-5, 'RelTol', 1e-5);
solinit = bvpinit(linspace(0, pi, 10), mat4init, lambda);
solinit = chebfun(mat4init, [0, pi]);

[sol, p] = bvp5c(mat4ode, mat4bc, solinit, lambda, opts);
sol1 = feval(sol, 1);
pass(2) = norm(sol1 - sol1_true, inf) < tol;
pass(3) = abs(p - p_true) < tol;


end

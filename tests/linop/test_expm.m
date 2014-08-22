function pass = test_expm
% TAD, 23 Jan 2014

tol = 1e-9; 
d = [-pi pi];
x = chebfun('x',d);
%%
[Z, I, D, C] = linop.primitiveOperators(d);
[z, E, s] = linop.primitiveFunctionals(d);
A = linop( D^2 );
A = addConstraint(A, E(-pi), 0);
A = addConstraint(A, E(pi), 0);

%%
% smooth initial condition
u0 = sin(exp(x)).*(pi^2-x.^2);
t = 0.02;
u = expm(A,t,u0);
exact = -4.720369127510475;
err(1,1) = abs( u(pi/2) - exact); 
err(1,2) = abs( u(-pi) );
err(1,3) = abs( u(pi) );

%%
% piecewise IC
tol = 1e-9; 
d = [-pi pi];
x = chebfun('x',d);
[Z, I, D, C] = linop.primitiveOperators(d);
[z, E, s] = linop.primitiveFunctionals(d);
A = linop( D^2 );
A = addConstraint(A, E(-pi), 0);
A = addConstraint(A, E(pi), 0);

u0 = chebfun(@(x) -abs(x)/pi+1, [-pi 0 pi]);
u = expm(A,0.01,u0);

exact = 0.95545945604534127;  % mathematica
err(2,1) = abs( u(.1) - exact);
err(2,2) = abs( u(-pi) );
err(2,3) = abs( u(pi) );

%%
% chebcolloc1
% smooth

prefs = cheboppref;
prefs.discretization = @chebcolloc1;
x = chebfun(@(x) x,d,'chebkind',1);
u0 = sin(exp(x)).*(pi^2-x.^2);
u = expm(A,0.02,u0,prefs);

exact = -4.720369127510475;
err(3,1) = abs( u(pi/2) - exact);
err(3,2) = abs( u(-pi) );
err(3,3) = abs( u(pi) );

%%
% piecewise IC
prefs.discretization = @chebcolloc1;
u0 = chebfun(@(x) -abs(x)/pi+1, [-pi 0 pi], 'chebkind', 1);
u = expm(A,0.01,u0,prefs);

exact = 0.95545945604534127;  % mathematica
err(4,1) = abs( u(.1) - exact);
err(4,2) = abs( u(-pi) );
err(4,3) = abs( u(pi) );

%%
% chebcolloc1
% smooth

prefs = cheboppref;
prefs.discretization = @ultraS;
x = chebfun(@(x) x,d);
u0 = sin(exp(x)).*(pi^2-x.^2);
u = expm(A,0.02,u0,prefs);

exact = -4.720369127510475;
err(5,1) = abs( u(pi/2) - exact);
err(5,2) = abs( u(-pi) );
err(5,3) = abs( u(pi) );

%%
% piecewise IC
prefs.discretization = @ultraS;
u0 = chebfun(@(x) -abs(x)/pi+1, [-pi 0 pi]);
u = expm(A,0.01,u0,prefs);

exact = 0.95545945604534127;  % mathematica
err(6,1) = abs( u(.1) - exact);
err(6,2) = abs( u(-pi) );
err(6,3) = abs( u(pi) );

%%

u0 = exp(-55*x.^2);
A = linop( D^2 + D );
A = addConstraint(A, E(-pi), 0);
A = addConstraint(A, E(pi), 0);
v = expm(A, 0, u0);
err(7, 1) = norm(u0 - v, inf);
err(7, 2) = length(u0) ~= length(v);

%%

pass = ( err < tol );

end

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
u0 = chebfun({[0;1] [1;0]},[-pi 0 pi]);
u = expm(A,0.01,u0);

exact = 0.95545945604534127;  % mathematica
err(2,1) = abs( u(.1) - exact);
err(2,2) = abs( u(-pi) );
err(2,3) = abs( u(pi) );

%%
% colloc1
% smooth

prefs = cheboppref;
prefs.discretization = @colloc1;
x = chebfun(@(x) x,d,'chebkind',1);
u0 = sin(exp(x)).*(pi^2-x.^2);
u = expm(A,0.02,u0,prefs);

exact = -4.720369127510475;
err(3,1) = abs( u(pi/2) - exact);
err(3,2) = abs( u(-pi) );
err(3,3) = abs( u(pi) );


%%
% piecewise IC
prefs.discretization = @colloc1;
x = chebfun(@(x) x,d,'chebkind',1);
u0 = min( 1-x/pi, 1+x/pi );
u = expm(A,0.01,u0,prefs);

exact = 0.95545945604534127;  % mathematica
err(4,1) = abs( u(.1) - exact);
err(4,2) = abs( u(-pi) );
err(4,3) = abs( u(pi) );


%%
pass = ( err < tol );

end

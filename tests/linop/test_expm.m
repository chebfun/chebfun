function pass = test_expm
% TAD, 23 Jan 2014

tol = 1e-9; 
d = [-pi pi];
x = chebfun('x',d);
pref = cheboppref();
pref.discretization = @chebcolloc2;
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
pref.discretization = @chebcolloc2;
u = expm(A,t,u0,pref);
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
u = expm(A,0.01,u0, pref);

exact = 0.95545945604534127;  % mathematica
err(2,1) = abs( u(.1) - exact);
err(2,2) = abs( u(-pi) );
err(2,3) = abs( u(pi) );

%%
% chebcolloc1
% smooth

pref.discretization = @chebcolloc1;
x = chebfun(@(x) x,d,'chebkind',1);
u0 = sin(exp(x)).*(pi^2-x.^2);
u = expm(A,0.02,u0,pref);

exact = -4.720369127510475;
err(3,1) = abs( u(pi/2) - exact);
err(3,2) = abs( u(-pi) );
err(3,3) = abs( u(pi) );

%%
% piecewise IC
pref.discretization = @chebcolloc1;
u0 = chebfun(@(x) -abs(x)/pi+1, [-pi 0 pi], 'chebkind', 1);
u = expm(A,0.01,u0,pref);

exact = 0.95545945604534127;  % mathematica
err(4,1) = abs( u(.1) - exact);
err(4,2) = abs( u(-pi) );
err(4,3) = abs( u(pi) );

%%
% chebcolloc1
% smooth

pref.discretization = @ultraS;
x = chebfun(@(x) x,d);
u0 = sin(exp(x)).*(pi^2-x.^2);
u = expm(A,0.02,u0,pref);

exact = -4.720369127510475;
err(5,1) = abs( u(pi/2) - exact);
err(5,2) = abs( u(-pi) );
err(5,3) = abs( u(pi) );

%%
% piecewise IC

u0 = chebfun(@(x) -abs(x)/pi+1, [-pi 0 pi]);
u = expm(A,0.01,u0,pref);

exact = 0.95545945604534127;  % mathematica
err(6,1) = abs( u(.1) - exact);
err(6,2) = abs( u(-pi) );
err(6,3) = abs( u(pi) );

%%

u0 = exp(-55*x.^2);
A = linop( D^2 + D );
A = addConstraint(A, E(-pi), 0);
A = addConstraint(A, E(pi), 0);
v = expm(A, 0, u0, pref);
err(7, 1) = norm(u0 - v, inf);
err(7, 2) = length(u0) ~= length(v);

%%

d = [-1 1];
[Z, I, D, C] = linop.primitiveOperators(d);
[z, E, s] = linop.primitiveFunctionals(d);
x = chebfun('x',d);
U0 = [ sin(2*pi*x).^2; 1 + 0*x ];
A = linop( [ D^2 Z; Z D^2 ] );
A = addConstraint(A, [ E(-1) z ], 0);
A = addConstraint(A, [ E(1) z ], 0);
A = addConstraint(A, [ z E(-1) ], 0);
A = addConstraint(A, [ z E(1) ], 0);
V = expm(A, 0, U0, pref);
err(8, 1) = norm(U0 - V, inf);
err(8, 2) = all(cellfun(@length,{U0{1:end}}) ~= ...
                cellfun(@length,{V{1:end}}));

%%

pass = ( err < tol );

end

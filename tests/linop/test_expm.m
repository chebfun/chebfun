function pass = test_expm
% TAD, 23 Jan 2014

tol = 1e-9; 
d = [-pi pi];
x = chebfun('x',d);

%%
[Z,I,D,C] = linop.primitiveOperators(d);
A = linop( D^2 );
A = addlbc(A,1);
A = addrbc(A,0);

%%
% smooth initial condition
u0 = sin(exp(x)).*(pi^2-x.^2) + (1-x/pi)/2;
t = 0.02;
u = expm(A,t,u0);
u = u{1};
exact = -4.47036912753658;  % mathematica
err(1) = abs( u(pi/2) - exact); 


%%
t = 0.005;
u0 = chebfun({[1;1] [1;0]},[-pi 0 pi]);
u = expm(A,t,u0);
u = u{1};

%%
exact = 0.99972973254164165445;  % mathematica
err(2) = abs( u(-.2) - exact);

%%
pass = ( err < tol );

end

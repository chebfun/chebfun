% A simple example BVP.
ccc

% The 'CHEBOP'
N.domain = [-1, 1];
N.op = @(x, u) diff(u, 2) + (x).*diff(u) + u;
N.lbc = @(u) u - 1;
N.bc = @(x, u) feval(u, 0.5);

%%
% Initialise an ADCHEBFUN:
x = chebfun('x', [-1 1]);
u = adchebfun(@(u) 0*u, [-1 1]);

%%
% Evaluate the operators to get a linearisation:
v = N.op(x, u);
A = v.jacobian;
Ares = v.func;
v = feval(N.lbc(u), -1);
B1 = v.jacobian;
B1res = v.func;
v = N.bc(x, u);
B2 = v.jacobian;
B2res = v.func;

%%

% Assemble to a CHEBMATRIX:
A = chebmatrix({A});
C = linopConstraint(chebmatrix({B1 ; B2}), -[B1res ; B2res]);
L = linop( A, C )


%%
% Make a CHEBMATRIX RHS:
rhs = chebmatrix( { x } );

%%
% Solve:
u = L\rhs;
u = u.blocks{1}
%%
% Display:

plot(u), shg
feval(N.lbc(u), N.domain(1))
N.bc(x, u)


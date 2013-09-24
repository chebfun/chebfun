% A simple example BVP.

% The 'CHEBOP':
N.op = @(x, u) diff(u, 2) + abs(x).*diff(u) + u;
N.lbc = @(u) u - 1;
N.bc = @(x, u) feval(u, 0.5);

%%
% Initialise an ADCHEBFUN:
x = chebfun('x', [-1 1]);
u = adchebfun(@(u) u, [-1 1]);

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
L = chebmatrix( { A } );
L.constraints(1).op = B1;
L.constraints(1).value = -B1res;
L.constraints(2).op = B2;
L.constraints(2).value = -B2res;

%%
% Make a CHEBMATRIX RHS:
rhs = chebmatrix( { x } );

%%
% Solve:
u = L\rhs;

%%
% Display:
u = u.blocks{1}
plot(u), shg
u(-1)-1
sum(u)
u(0)


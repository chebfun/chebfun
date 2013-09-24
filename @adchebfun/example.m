% A simple example BVP.

% The 'CHEBOP':
N.op = @(u) 0.1*diff(u, 2) - diff(u) + u;
N.lbc = @(u) u;
N.bc = @(u) sum(u);

% Initialise an ADCHEBFUN:
u = adchebfun(@(u) u);

% Evaluate the operators to get a linearisation:
Nu = N.op(u);
A = Nu.jacobian;
B1 = get(feval(N.lbc(u), -1), 'jacobian');
B2 = get(N.bc(u), 'jacobian');

% Assemble to a CHEBMATRIX:
L = chebmatrix( { A } );
L.constraints(1).op = B1;
L.constraints(1).value = 0;
L.constraints(2).op = B2;
L.constraints(2).value = 0;

% Make a CHEBMATRIX RHS:
rhs = chebmatrix( { chebfun('x') } );

% Solve:
u = L\rhs;

% Display:
u = u.blocks{1};
plot(u), shg
u(-1)
sum(u)





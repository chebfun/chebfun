clc
x = chebfun('x');
N = (1+linop.diag(x))*linop.diff([-1, 1], 2) + linop.diag(sin(x))*linop.diff + linop.diag(x)*linop.eye([-1, 1]); 
% N = linop.diff([-1, 1], 2) + linop.eye([-1,1]); 
B = chebmatrix({linop.feval(-1, [-1 1])});
B2 = chebmatrix({linop.feval(1, [-1 1])});

A = chebmatrix({N});
L = linop(A); 
bc = linopConstraint();
bc = append(bc, B, 1);
bc = append(bc, B2, 1);
L.constraint = bc;

% discretize(L.operator, 5)
% discretize(L.constraint.operator, 5)

f = chebfun(0);

tic
u = linsolve(L, f, @blockColloc2); 
toc
u = u{1};
u(0)
% u(.5)
plot(u), shg

tic
u = linsolve(L, f, @blockUS); 
toc
u = u{1};
u(0)
% u(.5)
hold on
plot(u,'r'), shg
hold off

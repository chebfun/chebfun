clc
x = chebfun('x');
N = (2+linop.diag(x.^2))*linop.diff([-1, 1], 2) + linop.diag(sin(100*x))*linop.diff + (1+linop.diag(x))*linop.eye([-1, 1]); 
% N = .0001*linop.diff([-1, 1], 2) + linop.diag(sin(x))*linop.diff + linop.diag(x)*linop.eye([-1, 1]); 
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

f = chebfun(@(x) cos(x));

tic
u = linsolve(L, f, @blockColloc2); 
toc
feval(u, 0)
plot(u{1}), shg

tic
v = linsolve(L, f, @blockUS); 
toc
feval(v, 0)
hold on
plot(v{1},'r'), shg
hold off

length(u{1})
length(v{1})
norm(u{1} - v{1})

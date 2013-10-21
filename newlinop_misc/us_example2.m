clc

dom = [-1, .5, 1];
x = chebfun('x', dom);
N = linop.diff(dom, 2) + linop.diag(sign(x))*linop.eye;

B = chebmatrix({linop.feval(-1, dom)});
B2 = chebmatrix({linop.feval(1, dom)});

A = chebmatrix({N});
L = linop(A); 
bc = linopConstraint();
bc = append(bc, B, 1);
bc = append(bc, B2, 1);
L.constraint = bc;

f = chebfun(@(x) cos(x), dom);

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
plot(v{1},'--r'), shg
hold off


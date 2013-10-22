clc

dom = [-1, 1];
x = chebfun('x', dom);
N = linop.diff(dom, 1) + linop.eye;
M = linop.diff(dom, 1) - linop.eye;
D = linop.diff(dom, 1) ;
Z = linop.zeros(dom);
I = linop.eye(dom);
A = chebmatrix({N D ; -D M})
discretize(A, 5)

z = linop.feval(dom(1), dom)*Z
B = chebmatrix({linop.feval(dom(1), dom), z});
B2 = chebmatrix({z, linop.feval(dom(2), dom)});

L = linop(A);
bc = linopConstraint();
bc = append(bc, B, 1);
bc = append(bc, B2, 1);
L.constraint = bc;

f = chebfun(@(x) 0*cos(10*x), dom);
f = [ f  ;f ];

tic
u = linsolve(L, f, @blockColloc2);
toc
feval(u, 0)
plot(u{1},'-b'); hold on
plot(u{2}, '.-b'), shg

tic
v = linsolve(L, f, @blockUS);
toc
feval(v, 0)

plot(v{1},'-r');
plot(v{2}, '.-r'), shg
hold off

norm(u{1} - v{1})
norm(u{2} - v{2})
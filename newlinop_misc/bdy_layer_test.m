ccc

I = linop.eye;
D = linop.diff;
D2 = linop.diff([-1 1], 2);
El = linop.feval(-1, [-1 1]);
Er = linop.feval(1, [-1 1]);
x = chebfun('x');

A = .01*D2 + D - I;
L = linop(A);
L = addbc(L, El, 1);
L = addbc(L, Er, 1);
L.discretization = @ultraS;

f = chebfun(0, [-1 -.5 1]);
u = L\f;


plot(u)
figure
chebpolyplot(u);

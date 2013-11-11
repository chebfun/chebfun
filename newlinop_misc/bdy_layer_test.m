ccc

I = linop.eye;
D = linop.diff;
D2 = linop.diff([-1 1], 2);
El = linop.feval(-1, [-1 1]);
Er = linop.feval(1, [-1 1]);

A = .001*D2 + D - I;
L = linop(A);
L.discretization = @ultraS;
L = addbc(L, El, 1);
L = addbc(L, Er, 1);

u = L\chebfun(0, [-1 0 1]);


u = u{1}
plot(u)
figure
chebpolyplot(u);
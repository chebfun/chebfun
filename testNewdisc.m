clear classes

%%
I = linop.eye();
D = linop.diff;
d = [-1 1];

%%
disc = colloc2(I);
disc.dimension = 6;
disc.domain = d;
discretize(disc)

%%
disc = colloc2(D^2+I);
disc.dimension = 6;
disc.domain = d;
discretize(disc)

%%
A = [ I D^2; -D D^2+I ];
disc = colloc2(A);
disc.dimension = [5 7];
disc.domain = [-1 0 1];
discretize(disc)

%%
L = linop(A);
x = chebfun('x');
f = [ x; sin(x) ];
El = linop.feval(-1,[-1 1]);
Er = linop.feval(1,[-1 1]);
z = linop.zero;
L = addbc(L,[El z],1);
L = addbc(L,[z Er],0);
L = addbc(L,[z El],-1);

d = colloc2(L);
d.domain = [-1 0 1];
d.dimension = [5 7];
d = deriveContinuity(d);

M = matrix(d)
rhs(d,f)


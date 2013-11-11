clear classes
clc

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
disc = colloc2(D^2);
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

%%
d = colloc2(L);
d.domain = [-1 0 1];
d.dimension = [5 7];
d = deriveContinuity(d);

M = matrix(d);
b = rhs(d,f)
size(M)




%%

ccc

%%
I = linop.eye();
D = linop.diff;
d = [-1 1];

%%
disc = ultraS(I);
disc.dimension = 6;
disc.domain = d;
full(discretize(disc))

%%
disc = ultraS(D^2+I);
disc.dimension = 6;
disc.domain = d;
full(discretize(disc))

%%
A = [ I D^2; -D D^2+I ];
disc = ultraS(A);
disc.dimension = [5 7];
disc.domain = [-1 0 1];
discretize(disc)

%%
clc
L = linop(A);
El = linop.feval(-1, [-1 1]);
Er = linop.feval(1, [-1 1]);
z = linop.zero;
Z = linop.zeros;
S = linop.sum([-1,1]);
L = addbc(L, [Er z], 0);
% L = addbc(L, [El z], 0);
L = addbc(L, [z El], 0);
L = addbc(L, [z S], 0);
L.discretization = @ultraS;
disc = ultraS(L);
% disc.dimension = [15 27];
% disc.domain = [-1 0 1];
% disc = deriveContinuity(disc);
% M = matrix(disc);
% spy(M), shg
% size(M)

f = chebfun(@sin);
u = L\[f ; f]


%%
ccc
x = chebfun('x');
I = linop.eye();
D = linop.diff;
El = linop.feval(-1, [-1 1]);
Er = linop.feval(1, [-1 1]);
z = linop.zero;
Z = linop.zeros;
S = linop.sum([-1,1]);

A = [ I D^2 exp(x); -D D^2+I x];
L = linop(A);

L = addbc(L, [Er z 0], 0);
% L = addbc(L, [El z], 0);
L = addbc(L, [z El 0], 0);
L = addbc(L, [z S 0], 0);
L.discretization = @ultraS;
disc = ultraS(L);
disc.dimension = [9];
% disc.domain = [-1 1];
disc.dimension = [15 27];
disc.domain = [-1 0 1];
disc = deriveContinuity(disc);
M = matrix(disc);
spy(M), shg
size(M)

%%

x = chebfun('x');
I = linop.eye();
D = linop.diff;
El = linop.feval(-1, [-1 1]);
Er = linop.feval(1, [-1 1]);
z = linop.zero;
Z = linop.zeros;
S = linop.sum([-1,1]);

A = [ I D^2 exp(x); -D D^2+I x ; S S 7];
L = linop(A);

L = addbc(L, [Er z 0], 0);
% L = addbc(L, [El z 1], 0);
L = addbc(L, [z El 0], 0);
L = addbc(L, [z S 0], 0);
L.discretization = @ultraS;
disc = ultraS(L);
disc.dimension = [100 100];
disc.domain = [-1 0 1];
disc = deriveContinuity(disc);
M = matrix(disc);
spy(M), shg
size(M)
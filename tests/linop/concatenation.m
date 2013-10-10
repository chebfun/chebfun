% Chebmatrix concatenation

clear classes
I = linop.eye([0 2]);
D = linop.diff([0 1 2]);
x = chebfun('x',[0 2]);
X = linop.diag(x);

%%
A = I + 2*X;
B = D^2;

C = [A,A];
size(matrix(C,5))

C = [A;A];
size(matrix(C,5))

%%
C = [A, B];
size(matrix(C,5))

C = [A; B];
size(matrix(C,5))


%%
C = [A, B];
E = [ C; 2*C ];
size(matrix(E,5))

C = [A, x];
size(matrix(C,5))

E = [C; C];
size(matrix(E,5))

%%
x = chebfun('x',[0 .5 2]);
size(matrix( [x,D], 5) )

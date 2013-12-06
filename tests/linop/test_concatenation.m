function pass = test_Concatnation

I = linop.eye([0 2]);
D = linop.diff([0 1 2]);
x = chebfun('x',[0 2]);
X = linop.mult(x);

%%
A = I + 2*X;
B = D^2;

C = [A,A];
C5 = cell2mat(discretize(C,5));
pass(1) = all(size(C5) == [5, 10]);

C = [A;A];
C5 = cell2mat(discretize(C,5));
pass(2) = all(size(C5) == [10, 5]);

%%
C = [A, B];
C5 = cell2mat(discretize(C,5));
pass(3) = all(size(C5) == [10, 20]);

C = [A; B];
C5 = cell2mat(discretize(C,5));
pass(4) = all(size(C5) == [20, 10]);

%%
C = [A, B];
E = [ C; 2*C ];
C5 = cell2mat(discretize(C,5));
pass(5) = all(size(C5) == [10, 20]);

C = [A, x];
C5 = cell2mat(discretize(C,5));
pass(6) = all(size(C5) == [5,  6]);

E = [C; C];
C5 = cell2mat(discretize(C,5));
pass(7) = all(size(C5) == [5,  6]);

%%
x = chebfun('x',[0 .5 2]);
C5 = cell2mat(discretize([x, D],5));
pass(8) = all(size(C5) == [15,  16]);


end
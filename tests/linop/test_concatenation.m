function pass = test_concatenation

I = operatorBlock.eye([0 2]);
D = operatorBlock.diff([0 1 2]);
x = chebfun('x',[0 2]);
X = operatorBlock.mult(x);

%%
A = I + 2*X;
B = D^2;

C = [A,A];
C5 = matrix(C,5);
pass(1) = all(size(C5) == [5, 10]);

C = [A;A];
C5 = matrix(C,5);
pass(2) = all(size(C5) == [10, 5]);

%%
C = [A, B];
C5 = matrix(C,[5 5]);
pass(3) = all(size(C5) == [10, 20]);

C = [A; B];
C5 = matrix(C,[5 5]);
pass(4) = all(size(C5) == [20, 10]);

%%
C = [A, B];
E = [ C; 2*C ];
C5 = matrix(C,[5 5]);
pass(5) = all(size(C5) == [10, 20]);

C = [A, x];
C5 = matrix(C,5);
pass(6) = all(size(C5) == [5,  6]);

E = [C; C];
C5 = (matrix(C,5));
pass(7) = all(size(C5) == [5,  6]);

%%
x = chebfun('x',[0 .5 2]);
C5 = (matrix([x, D],[5 5 5]));
pass(8) = all(size(C5) == [15,  16]);


end

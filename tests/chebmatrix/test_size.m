function pass = test_size(pref)
%TEST_SIZE     Test the SIZE method of the CHEBMATRIX class
%%
% Create some CHEBMATRIX objects and check their sizes:
x = chebfun(@(x) x);
D = operatorBlock.diff();
S = functionalBlock.sum();

ff = chebmatrix([x, x, x]);
A = [[D, D, D]; [S, S, S]];
B = [[D x]; [S 1]];

pass(1) = all(size(ff) == [1 3]);
pass(2) = all(size(ff') == [3 1]);
pass(3) = all(size(A) == [2 3]);
pass(4) = all(size(A') == [3 2]);
pass(5) = all(size(B) == [2 2]);
pass(6) = all(size(B) == [2 2]);
%%
end


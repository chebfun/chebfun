function pass = test_length(pref)
%TEST_SIZE     Test the LENGTH method of the CHEBMATRIX class
%%
% Create some CHEBMATRIX objects and check their lengths:
x = chebfun(@(x) x);
D = operatorBlock.diff();
S = functionalBlock.sum();

ff = chebmatrix([x, x, x]);
A = [[D, D, D]; [S, S, S]];
B = [[D x]; [S 1]];

pass(1) = all(length(ff) == 3);
pass(2) = all(length(ff') == 3);
pass(3) = all(length(A) == 3);
pass(4) = all(length(A') == 3);
pass(5) = all(length(B) == 2);
pass(6) = all(length(B) == 2);
%%
end


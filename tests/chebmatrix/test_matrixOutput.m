function pass = test_matrixOutput(~)
% TEST_MATRIXOUTPUT   Chebmatrix operations resulting only in doubles should
% return a normal matrix.

dom = [-3 1];
% A chebmatrix with 5 column functions
T = chebpoly(1:5, dom);
T = chebmatrix(T);
C = T'*T;
pass(1) = isnumeric(C) && all( size(C) == [5 5]);

% Two rows, 5 columns
TT = [T; T];
CC = TT'*TT;
pass(2) = isnumeric(CC) && all( size(CC) == [5 5]);

% Operators involved, output should not be numeric
[Z, I, D, C, M] = linop.primitiveOperators(dom);
A = [I D; Z C];
f = [chebpoly(1, dom); chebpoly(5, dom)];
Af = A*f;
pass(3) = isa(Af, 'chebmatrix') && ...
    all(cellfun(@isa, Af.blocks, {'chebfun'; 'chebfun'})) && ...
    all( size(Af) == [2 1]);
end
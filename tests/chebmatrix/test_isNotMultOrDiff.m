function pass = test_isNotMultOrDiff(~)
%TEST_ISNOTMULTORDIFF    Test that the property of operator not being
%diff/integration is dealt with correctly.

% Construct the primitive functionalBlocks
dom = [0, 2];
[Z, E, S, D] = linop.primitiveFunctionals(dom);
E = E(1);
D = D(chebfun(@sin, dom));

% Throw in some operatorBlocks for good measure 
[ZZ, I, DD, C, M] = linop.primitiveOperators(dom);
M = M(chebfun(@sin, dom));

% Concatenate to get chebmatrices, and check that they have the expected value
% of the property
A = [ZZ DD; C M];
pass(1) = all(all(A.isNotDiffOrInt == [1 0; 0 1]));
B = [I; S];
pass(2) = all(all(B.isNotDiffOrInt == [1; 0]));
C = 2*A;
pass(3) = all(all(C.isNotDiffOrInt == [1 0; 0 1]));
D = B + B;
pass(4) = all(all(D.isNotDiffOrInt == [1; 0]));
E = A*A;
pass(5) = all(all(E.isNotDiffOrInt == [0 0; 0 0]));

end
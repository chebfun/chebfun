function pass = test_isNotMultOrDiff(~)
%TEST_ISNOTMULTORDIFF    Test that the property of operator not being
%diff/integration is dealt with correctly.

% Construct the primitive functionalBlocks
dom = [0, 2];
[Z, E, S, D] = linop.primitiveFunctionals(dom);
E = E(1);
D = D(chebfun(@sin, dom));

% Throw in some operatorBlocks for good measure (tested later in combinations
% with functionalBlocks)
[ZZ, I, DD, C, M] = linop.primitiveOperators(dom);
M = M(chebfun(@sin, dom));

% Check that they have the expected value of the property
pass(1) = ( Z.isNotDiffOrInt == 1 );
pass(2) = ( E.isNotDiffOrInt == 1 );
pass(3) = ( S.isNotDiffOrInt == 0 );
pass(4) = ( D.isNotDiffOrInt == 1 );

% Operations on functionalBlocks
A = Z + I;
pass(5) = ( A.isNotDiffOrInt == 1 );
A = E + S;
pass(6) = ( A.isNotDiffOrInt == 0 );
A = 2*Z;
pass(7) = ( A.isNotDiffOrInt == 1 );
A = 2*S;
pass(8) = ( A.isNotDiffOrInt == 0 );
A = S*ZZ;
pass(9) = ( A.isNotDiffOrInt == 1 );
A = D*DD;
pass(10) = ( A.isNotDiffOrInt == 0 );
A = E*C;
pass(11) = ( A.isNotDiffOrInt == 0 );
end
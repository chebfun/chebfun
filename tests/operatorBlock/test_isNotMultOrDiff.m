function pass = test_isNotMultOrDiff(~)
%TEST_ISNOTMULTORDIFF    Test that the property of operator not being
%diff/integration is dealt with correctly.

% Construct the primitive operatorBlocks
dom = [0, 2];
[Z, I, D, C, M] = linop.primitiveOperators(dom);
M = M(chebfun(@sin, dom));

% Check that they have the expected value of the property
pass(1) = ( Z.isNotDiffOrInt == 1 );
pass(2) = ( I.isNotDiffOrInt == 1 );
pass(3) = ( D.isNotDiffOrInt == 0 );
pass(4) = ( C.isNotDiffOrInt == 0 );
pass(5) = ( M.isNotDiffOrInt == 1 );

% Operations on operatorBlocks
A = Z + I;
pass(6) = ( A.isNotDiffOrInt == 1 );
A = Z + D;
pass(7) = ( A.isNotDiffOrInt == 0 );
A = M*Z;
pass(8) = ( A.isNotDiffOrInt == 1 );
A = C*Z;
pass(9) = ( A.isNotDiffOrInt == 1 );
A = C + M;
pass(10) = ( A.isNotDiffOrInt == 0 );
A = 2*M;
pass(11) = ( A.isNotDiffOrInt == 1 );
A = 2*D;
pass(12) = ( A.isNotDiffOrInt == 0 );
end
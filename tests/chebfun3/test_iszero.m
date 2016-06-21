function pass = test_iszero()
% Test chebfun3/iszero.

f = chebfun3(0);
pass(1) = iszero(f);

f = chebfun3([]);
pass(2) = iszero(f);

f = chebfun3(2);
pass(3) = ~iszero(f);

f = chebfun3(@(x,y,z) x+y-z);
pass(4) = ~iszero(f);
end
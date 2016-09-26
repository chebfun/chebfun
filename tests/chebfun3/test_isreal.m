function pass = test_isreal()
% Test chebfun3/isreal.

f = chebfun3(@(x,y,z) x+y-z);
pass(1) = isreal(f);

f = chebfun3(@(x,y,z) 1i*x+y-z);
pass(2) = ~isreal(f);

f = real(chebfun3(@(x,y,z) 1i*x+y-z));
pass(3) = isreal(f);

end
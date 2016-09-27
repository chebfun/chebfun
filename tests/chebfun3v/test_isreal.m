function pass = test_isreal()
% Test chebfun3v/isreal.

f = chebfun3v(@(x,y,z) x, @(x,y,z) y - z);
pass(1) = isreal(f);

f = chebfun3v(@(x,y,z) 1i*x, @(x,y,z) y - z);
pass(2) = ~isreal(f);

f = real(chebfun3v(@(x,y,z) 1i*x, @(x,y,z) y - z));
pass(3) = isreal(f);

end
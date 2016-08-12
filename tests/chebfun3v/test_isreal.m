function pass = test_isreal()
% Test chebfun3v/isreal.

j = 1;

f = chebfun3v(@(x,y,z) x, @(x,y,z) y - z);
pass(j) = isreal(f);
j = j+1;

f = chebfun3v(@(x,y,z) 1i*x, @(x,y,z) y - z);
pass(j) = ~isreal(f);
j = j+1;

f = real(chebfun3v(@(x,y,z) 1i*x, @(x,y,z) y - z));
pass(j) = isreal(f);

end
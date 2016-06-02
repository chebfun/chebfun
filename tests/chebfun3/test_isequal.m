function pass = test_isequal()
% Test @chebfun3/isequal.m
% [reviewed by LNT 31.05.16]

f = chebfun3(@(x,y,z) cos(x.*y.*z));
g = chebfun3(@(x,y,z) sin(x + y.^2 - z.^3));

pass(1) = isequal(f, f);
pass(2) = ~isequal(f, g);

end
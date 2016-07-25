function pass = test_domainChck()
% Test chebfun3/domainCheck

% Default domain: 
ff = @(x,y,z) cos(x+y+z);
gg = @(x,y,z) sin(x.*y.*z);
f = chebfun3(ff);
g = chebfun3(gg);
pass(1) = domainCheck(f, g);

% A different domain: 
dom = [-1 2 -2 1 -3 0];
ff = @(x,y,z) cos(x+y+z);
gg = @(x,y,z) sin(x.*y.*z);
f = chebfun3(ff, dom);
g = chebfun3(gg, dom);
pass(2) = domainCheck(f, g);

end 
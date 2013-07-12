a = .5; b = 1.5;
op = @(x) sin(1-x)./((1-x).^a.*(1+x).^b);
exponents = [];
isSingEnd = [1,1];
singType = {'branch', 'branch'};
pref = [];
f = singfun( op, exponents, isSingEnd, singType, pref )

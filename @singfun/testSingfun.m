a = .5; b = .5;
f = singfun( @(x) sin(x)./((1-x).^a.*(1+x).^b), [1 1], {'branch', 'branch'} )
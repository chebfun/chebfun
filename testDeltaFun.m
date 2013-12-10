deltafun

%%
% shouldn't work
deltafun(1)

%%
deltafun(1,2)
%%
% shouldn't work
deltafun( rand(3,1), rand(1,4) )
%%
% shouldn't work
deltafun( rand(3,1), rand(3,1) )
%%
% should work
deltafun( rand(3,3), rand(3,1) )
deltafun( rand(3,3), rand(1,3) )
%%
% shouldn't work
deltafun( rand(3,1), rand(1,3) )
%%
% should work
deltafun( rand(1,3), rand(3,1) )
%%
% Test simplify
mag = [rand(3,3); zeros(3,3)]; loc = rand(1,3);
%mag(end,end) = 1e-12;
d = deltafun(mag, loc, chebfun(0)); d.delta.magnitude
d = simplify(d); d.delta.magnitude

%%
mag = [rand(3,3); zeros(3,3)]; loc = [-1 -1+4*eps -1+8*eps ]; 
mag(end,end) = 1e-12;
d = deltafun(mag, loc, chebfun(0)); d.delta.magnitude
d = simplify(d); d.delta.magnitude - sum(mag,2)


%%
% Dirac delta and derivatives
d = dirac(deltafun);
d = diff(d,5);
d.delta.magnitude
%%
% The dirac delta function and inner products
loc = 0;
mag = 1;
f = chebfun(0);
d = deltafun(mag, loc, f)
x = chebfun('x')
ip(d, x)
f = sin(x);
ip(diff(d), f )

%%
mag = rand(1,5);
loc = rand(1,5);
d = deltafun( mag, loc, chebfun(0));
ip( d, sin(2*pi*x))
%%
% Some more inner products
loc = rand(1, 5);
mag = rand(1, 5);
d = deltafun( mag, loc, chebfun(0))
ip(d, 1) - sum(mag)
%%
loc = rand(1, 4);
mag = rand(1, 4);
x = chebfun('x');
d = deltafun( mag, loc, 0*x)
ip(d, x) - sum(mag.*loc)

%%
n = 4
ip( diff(dirac(deltafun),n), x.^n ) - (-1)^n * factorial(n)

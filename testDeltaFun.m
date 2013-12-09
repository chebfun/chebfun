deltafun

%%
deltafun(1)

%%
deltafun(1,2)
%%
deltafun( rand(3,1), rand(1,4) )
%%
deltafun( rand(3,1), rand(3,1) )
deltafun( rand(3,1), rand(1,3) )
deltafun( rand(1,3), rand(3,1) )
%%
% The dirac delta function and inner products
loc = 0;
mag = 1;
f = chebfun(0);
d = deltafun(mag, loc, f)
x = chebfun('x')
ip(d, x)
f = cos(x);
ip(d, f )

%%
% Some more inner products
loc = rand(5,1);
mag = rand(5,1);
d = deltafun( mag, loc, chebfun(0))
ip(d, 1) - sum(mag)
%%
loc = rand(5,1);
mag = rand(5,1);
x = chebfun('x');
d = deltafun( mag, loc, 0*x)
ip(d, x) - sum(mag.*loc)

%%

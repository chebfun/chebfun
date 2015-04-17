


%% Construction: 
% A few tests for construction: 

f = @(x,y,z) x.^2 + y.^2 + z.^2;
g = spherefun( f )

f = @(x,y,z) exp(-cos(pi*(x+y+z)));
g = spherefun( f )

f = @(x,y,z) exp(x);
g = spherefun( f )

f = @(x,y,z) exp(y);
g = spherefun( f )

f = @(x,y,z) exp(z);
g = spherefun( f )

f = @(x,y,z) cos(x.*y);
g = spherefun( f )

f = @(x,y,z) sin(x.*y.*z);
g = spherefun( f )

f = @(x,y,z) sin(x+ y.*z);
g = spherefun( f )

f = @(x,y,z) sin(x+ y.*z) + 1;
g = spherefun( f )


%% Pointwise evaluation: 

lam = rand; th = rand; 

x = cos(lam).*cos(th);
y = sin(lam).*cos(th);
z = sin(th);

f = @(x,y,z) sin(x+ y.*z) + 1;
g = spherefun( f );
feval(g, lam, th) - f(x,y,z)









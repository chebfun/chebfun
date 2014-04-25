deltafun

%%
% shouldn't work
deltafun(1)

%%
deltafun(1,2)
%% 
deltafun(1,0)
%%
% shouldn't work
%deltafun( rand(3,1), rand(1,4) )
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
% constructor simplifies
d = deltafun([], mag, loc); d.impulses
d.deltaLoc

%%
mag = [rand(3,3); zeros(3,3)]; loc = [-1 -1+4*eps -1+8*eps ]; 
mag(end,end) = 1e-12;
d = deltafun([], mag, loc); d.impulses, d.deltaLoc
d = simplify(d); d.impulses - sum(mag,2)


%%
% Dirac delta and derivatives
d = dirac(deltafun);
d = diff(d,5);
d.impulses
%%
% The dirac delta function and inner products
loc = 0;
mag = 1;
T = chebfun(0);
d = deltafun(mag, loc, T)
x = chebfun('x')
ip(d, x)
T = sin(x);
ip(diff(d), T )

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

%%
x = chebfun('x');
n = 5;
mag = 5*rand(1,n);
loc = linspace(-1,1, n + 2);
loc = loc(2:end-1);
d = deltafun([], mag, loc);

plot(d)

%%
% The Taylor series operator
N = 7;
x = chebfun('x');
h = .1;
x0 = 0;
d = dirac(x-x0);
T = 0;
for i = 0:N
    T = T + (-1)^i*d.*h^i/factorial(i);
    d = diff(d);
end
T = T.funs{1};
f = exp(x) + sin(x.^2);
innerProduct(T, f )-f(x0+h)
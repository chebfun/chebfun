function pass = test_atan2(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

a = -2.25*pi;
b = 2.25*pi;

%% Scalar-valued, tan(f, g):
x = chebfun(@(x) x, [a, b], pref);
f = .5+sin(x).*exp(-.1*x.^2);
g = cos(x).*(1+x.^2);
h = atan2(f, g);
tol = 10*epslevel(h);

xx = linspace(.99*a, .99*b, 100);
ff = feval(f, xx);
gg = feval(g, xx);
hh = atan2(ff, gg);
pass(1) = norm(feval(h, xx) - hh, inf) < 3*tol;

ends = h.domain([1 end]);
fi = feval(f, ends);
gi = feval(g, ends);
hi = atan2(fi, gi).';
pass(2) = norm(hi - h.pointValues([1 end])) < 2*tol;

%% Scalar-valued, tan(g, f):
h = atan2(g, f);
tol = 10*epslevel(h).*vscale(h);
hh = atan2(gg, ff);
err = norm(feval(h, xx) - hh, inf);
pass(3) = err < tol;

ends = h.domain([1 end]);
hi = atan2(feval(g, ends), feval(f, ends)).';
pass(4) = norm(hi - h.pointValues([1 end])) < tol;

%% y has a zero FUN:
x = chebfun('x');
f = 0*x;
g = sin(4*x);
h = atan2(f, g);
tol = 10*epslevel(h);

xx = linspace(-.99, .99, 100);
ff = feval(f, xx);
gg = feval(g, xx);
hh = atan2(ff, gg);
pass(5) = norm(feval(h, xx) - hh, inf) <= 10*vscale(h)*epslevel(h);

%% x has a zero FUN:
h = atan2(g, f);
tol = 10*epslevel(h);

hh = atan2(gg, ff);
pass(6) = norm(feval(h, xx) - hh, inf) <= 10*vscale(h)*epslevel(h);

end
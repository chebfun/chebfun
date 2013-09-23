function pass = test_atan2(pref)

if ( nargin == 0 )
    pref = chebfun.pref();
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
pass(1) = norm(feval(h, xx) - hh, inf) < tol;

ends = h.domain;
fi = feval(f, ends);
gi = feval(g, ends);
hi = atan2(fi, gi).';
pass(2) = norm(hi - h.impulses(:,:,1)) < tol;

%% Scalar-valued, tan(g, f):
h = atan2(g, f);
tol = 10*epslevel(h);
hh = atan2(gg, ff);
pass(3) = norm(feval(h, xx) - hh, inf) < tol;

ends = h.domain;
hi = atan2(feval(g, ends), feval(f, ends)).';
pass(4) = norm(hi - h.impulses(:,:,1)) < tol;

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
pass(5) = norm(feval(h, xx) - hh, inf);
ends = h.domain;
fi = feval(f, ends);
fi(abs(fi) < tol) = 0;
gi = feval(g, ends);
gi(abs(gi) < tol) = 0;
hi = atan2(fi, gi).';
pass(6) = norm(hi - h.impulses(:,:,1)) < epslevel(h);

%% x has a zero FUN:
h = atan2(g, f);
tol = 10*epslevel(h);

hh = atan2(gg, ff);
pass(7) = norm(feval(h, xx) - hh, inf);

ends = h.domain;
fi = feval(f, ends);
fi(abs(fi) < tol) = 0;
gi = feval(g, ends);
gi(abs(gi) < tol) = 0;
hi = atan2(gi, fi).';
pass(8) = norm(hi - h.impulses(:,:,1)) < tol;

end


function pass = test_hypot(pref)

if ( nargin == 0 ) 
    pref = chebfunpref();
end

% TODO: This test is meaningless

% Test the example given in the file.
f = chebfun(@(x) 3*[1e300*x 1e-300*x], pref);
g = chebfun(@(x) 4*[1e300*x 1e-300*x], pref);
h = hypot(f, g);

seedRNG(6178);
xx = 2 * rand(100, 1) - 1;
ff = 3*[1e300*xx 1e-300*xx];
gg = 4*[1e300*xx 1e-300*xx];
hh = hypot(ff, gg);
pass(1) = norm(feval(h, xx) - hh, inf)/vscale(h) < 10*eps;

% Test a function with breakpoints.
pref.splitting = 1;
base_op = @(x) sign(x - 0.1).*abs(x + 0.2).*sin(3*x)*(pi/16) + pi/8;
f = chebfun(base_op, [-1 -0.2 0.1 1], pref);
g = hypot(f, f);
g_op = @(x) hypot(base_op(x), base_op(x));
pass(2) = norm(feval(g, xx) - g_op(xx), inf) < 1e2*vscale(g)*eps;
    

end

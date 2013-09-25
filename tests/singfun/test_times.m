% Test file for singfun/times.m

function pass = test_times(pref)

% Get preferences.
if ( nargin < 1 )
    pref = singfun.pref;
end

% Generate a few random points to use as test values.
seedRNG(6178);
d = 14;
x = 2*(1-10^(-d)) * rand(500, 1) - (1-10^(-d));

% The order of the exponents:
a = 0.64;
b = -0.64;
c = 1.28;
d = -1.28;
p = -0.2;
q = -0.3;

% Pre-allocate pass matrix
pass = zeros(1, 7);

% fractional pole at the right endpoint
f = singfun(@(x) (1+x).^p, [p 0], {'sing', 'none'}, [], [], pref);
g = singfun(@(x) (1+x).^q, [q 0], {'sing', 'none'}, [], [], pref);
h = f.*g;
vals_h = feval(h, x);
h_exact = @(x) 1./sqrt(1+x);
vals_exact = feval(h_exact, x);
err = vals_h - vals_exact;
pass(1) = all( abs(err) < max(get(f, 'epslevel'), get(g, 'epslevel'))*abs(vals_exact) );

% root at the left endpoint
f = singfun(@(x) (1+x).^c.*sin(x), [c 0], {'root', 'none'}, [], [], pref);
g = singfun(@(x) (1+x).^(2*c), [2*c 0], {'root', 'none'}, [], [], pref);
h = f.*g;
vals_h = feval(h, x);
h_exact = @(x) (1+x).^(3*c).*sin(x);
vals_exact = feval(h_exact, x);
err = vals_h - vals_exact;
pass(2) = all( abs(err) < 1e1*max(get(f, 'epslevel'), get(g, 'epslevel'))*abs(vals_exact) );

% fractional root at the right endpoint
f = singfun(@(x) (1-x).^c.*cos(x), [0 c], {'none', 'root'}, [], [], pref);
g = singfun(@(x) (1-x).^a.*cos(x), [0 a], {'none', 'root'}, [], [], pref);
h = f.*g;
vals_h = feval(h, x);
h_exact = @(x) (1-x).^(a+c).*(cos(x).^2);
vals_exact = feval(h_exact, x);
err = vals_h - vals_exact;
pass(3) = all( abs(err) < 1e1*max(get(f, 'epslevel'), get(g, 'epslevel'))*abs(vals_exact) );

% fractional pole at the right endpoint
f = singfun(@(x) (1-x).^b.*(x.^5), [0 b], {'none', 'sing'}, [], [], pref);
g = singfun(@(x) exp(x).*sin(5*x), [0 0], {'none', 'none'}, [], [], pref);
h = f.*g;
vals_h = feval(h, x);
h_exact = @(x) (1-x).^b.*(x.^5).*exp(x).*sin(5*x);
vals_exact = feval(h_exact, x);
err = vals_h - vals_exact;
pass(4) = all( abs(err) < 1e2*max(get(f, 'epslevel'), get(g, 'epslevel')) );

% a combination of fractional pole and fractional root
f = singfun(@(x) (1+x).^b.*sin(x), [b 0], {'sing', 'none'}, [], [], pref);
g = singfun(@(x) sin(2*x).*(1-x).^c, [0 c], {'none', 'root'}, [], [], pref);
h = f.*g;
vals_h = feval(h, x);
h_exact = @(x) (1+x).^b.*sin(x).*sin(2*x).*(1-x).^c;
vals_exact = feval(h_exact, x);
err = vals_h - vals_exact;
pass(5) = all( abs(err) < 1e2*max(get(f, 'epslevel'), get(g, 'epslevel')) );

% poles at different endpoints
f = singfun(@(x) sin(x).*(1-x.^2).^b, [b b], {'sing', 'sing'}, [], [], pref);
g = singfun(@(x) cos(x).^3.*(1+x).^p, [p 0], {'sing', 'none'}, [], [], pref);
h = f.*g;
vals_h = feval(h, x);
h_exact = @(x) sin(x).*(1-x).^b.*cos(x).^3.*(1+x).^(b+p);
vals_exact = feval(h_exact, x);
err = vals_h - vals_exact;
pass(6) = all( abs(err) < 1e1*max(get(f, 'epslevel'), get(g, 'epslevel'))*abs(vals_exact) );

% Check the trivial case with both vanishing alpha and beta.
f = singfun(@(x) exp(x).*x.^3.*sin(2*x), [0 0], {'none', 'none'}, [], [], pref);
g = singfun(@(x) exp(1-x).^(3/2), [0 0], {'none', 'none'}, [], [], pref);
h = f.*g;
vals_h = feval(h, x);
h_exact = @(x) exp(x).*x.^3.*sin(2*x).*exp(1-x).^(3/2);
vals_exact = feval(h_exact, x);
err = vals_h - vals_exact;
pass(7) = all( abs(err) < 1e2*max(get(f, 'epslevel'), get(g, 'epslevel')) );

end
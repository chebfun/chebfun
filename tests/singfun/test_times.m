% Test file for singfun/times.m

function pass = test_times(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
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

% Check operation in the case of empty arguments.
f = singfun();
data.exponents = [-1, 0];
g = singfun(@(x) 1./(1+x), data);
pass(1) = (isempty(f.*f) && isempty(f.*g) && isempty(g.*f));

% Multiplication of smooth SINGFUNs should not return a SINGFUN
f = singfun(@(x) sin(x));
g = singfun(@(x) cos(x));
pass(2) = ~isa(f.*g, 'singfun');

% SMOOTHFUN .* SINGFUN
f = smoothfun.constructor(@(x) sin(x));
g = singfun(@(x) cos(x));
pass(3) = isa(f.*g, 'smoothfun') && isa(g.*f, 'smoothfun');

% fractional pole at the right endpoint
data.exponents = [p 0];
data.singType = {'sing', 'none'};
f = singfun(@(x) (1+x).^p, data, pref);
data.exponents = [q 0];
data.singType = {'sing', 'none'};
g = singfun(@(x) (1+x).^q, data, pref);
h = f.*g;
vals_h = feval(h, x);
h_exact = @(x) 1./sqrt(1+x);
vals_exact = feval(h_exact, x);
err = vals_h - vals_exact;
pass(4) = all( abs(err) < 10*eps* ...
    abs(vals_exact) );

% root at the left endpoint
data.exponents = [c 0];
data.singType = {'root', 'none'};
f = singfun(@(x) (1+x).^c.*sin(x), data, pref);
data.exponents = [2*c 0];
data.singType = {'root', 'none'};
g = singfun(@(x) (1+x).^(2*c), data, pref);
h = f.*g;
vals_h = feval(h, x);
h_exact = @(x) (1+x).^(3*c).*sin(x);
vals_exact = feval(h_exact, x);
err = vals_h - vals_exact;
bound = 1e1*eps*norm(vals_exact, inf);
pass(5) = ( norm(err, inf) < bound );

% fractional root at the right endpoint
data.exponents = [0 c];
data.singType = {'none', 'root'};
f = singfun(@(x) (1-x).^c.*cos(x), data, pref);
data.exponents = [0 a];
data.singType = {'none', 'root'};
g = singfun(@(x) (1-x).^a.*cos(x), data, pref);
h = f.*g;
vals_h = feval(h, x);
h_exact = @(x) (1-x).^(a+c).*(cos(x).^2);
vals_exact = feval(h_exact, x);
err = vals_h - vals_exact;
bound = 1e1*eps*norm(vals_exact, inf);
pass(6) = ( norm(err, inf) < bound );

% fractional pole at the right endpoint
data.exponents = [0 b];
data.singType = {'none', 'sing'};
f = singfun(@(x) (1-x).^b.*(x.^5), data, pref);
data.exponents = [0 0];
data.singType = {'none', 'none'};
g = singfun(@(x) exp(x).*sin(5*x), data, pref);
h = f.*g;
vals_h = feval(h, x);
h_exact = @(x) (1-x).^b.*(x.^5).*exp(x).*sin(5*x);
vals_exact = feval(h_exact, x);
err = norm(vals_h - vals_exact, inf);
tol = 10*norm(vals_h, inf)*eps;
pass(7) = err < tol;

% a combination of fractional pole and fractional root
data.exponents = [b 0];
data.singType = {'sing', 'none'};
f = singfun(@(x) (1+x).^b.*sin(x), data, pref);
data.exponents = [0 c];
data.singType = {'none', 'root'};
g = singfun(@(x) sin(2*x).*(1-x).^c, data, pref);
h = f.*g;
vals_h = feval(h, x);
h_exact = @(x) (1+x).^b.*sin(x).*sin(2*x).*(1-x).^c;
vals_exact = feval(h_exact, x);
err = vals_h - vals_exact;
pass(8) = all( abs(err) < 1e4*eps );
    
    
% poles at different endpoints
data.exponents = [b b];
data.singType = {'sing', 'sing'};
f = singfun(@(x) sin(x).*(1-x.^2).^b, data, pref);
data.exponents = [p 0];
data.singType = {'sing', 'none'};
g = singfun(@(x) cos(x).^3.*(1+x).^p, data, pref);
h = f.*g;
vals_h = feval(h, x);
h_exact = @(x) sin(x).*(1-x).^b.*cos(x).^3.*(1+x).^(b+p);
vals_exact = feval(h_exact, x);
err = norm((vals_h - vals_exact)./(vals_exact), inf);
tol = 1e3*eps;

pass(9) = all( err < tol );

% Check the trivial case with both vanishing alpha and beta.
data.exponents = [0 0];
data.singType = {'none', 'none'};
f = singfun(@(x) exp(x).*x.^3.*sin(2*x), data, pref);
data.exponents = [0 0];
data.singType = {'none', 'none'};
g = singfun(@(x) exp(1-x).^(3/2), data, pref);
h = f.*g;
vals_h = feval(h, x);
h_exact = @(x) exp(x).*x.^3.*sin(2*x).*exp(1-x).^(3/2);
vals_exact = feval(h_exact, x);
err = norm(vals_h - vals_exact, inf);
tol = 1e3*norm(vals_h, inf)*eps;
pass(10) = err < tol;

end

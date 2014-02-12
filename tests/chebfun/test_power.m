function pass = test_power(pref)

if ( nargin == 0 ) 
    pref = chebpref();
end

% [TODO]: So far, POWER only supports easy cases (i.e., positive integers).
% Tests will have to expanded once this functionality is available.


%% Scalar-valued
f = chebfun(@(x) sin(x));
g = f.^0;
pass(1) = normest(g - 1) < 10*epslevel(f);

g = f.^1;
pass(2) = normest(g - f) < 10*epslevel(f);

g = f.^2;
h = chebfun(@(x) sin(x).^2);
pass(3) = normest(g - h) < 10*epslevel(h);

g = f.^3;
h = chebfun(@(x) sin(x).^3);
pass(4) = normest(g - h) < 10*epslevel(h);

%% Array-valued
f = chebfun(@(x) [sin(x), cos(x), 1i*exp(x)]);
fq = quasimatrix(f);

g = f.^0;
pass(5) = min(size(g)) == 3 && normest(g - 1) < epslevel(f);
gq = fq.^0;
pass(6) = normest(gq - g) < epslevel(g)*vscale(g);

g = f.^1;
pass(7) = min(size(g)) == 3 && normest(g - f) < epslevel(f);
gq = fq.^1;
pass(8) = normest(gq - g) < epslevel(g)*vscale(g);

g = f.^2;
h = chebfun(@(x) [sin(x).^2, cos(x).^2, -exp(2*x)]);
pass(9) = min(size(g)) == 3 && normest(g - h) < 10*vscale(h)*epslevel(h);
gq = fq.^2;
pass(10) = normest(gq - g) < epslevel(g);

g = f.^3;
h = chebfun(@(x) [sin(x).^3, cos(x).^3, -1i*exp(3*x)]);
pass(11) = min(size(g)) == 3 && normest(g - h) < 10*vscale(h)*epslevel(h);
gq = fq.^3;
pass(12) = normest(gq - g) < epslevel(g)*vscale(g);

%% constant .^ CHEBFUN
f = chebfun(@(x) sin(x));
g = 1.^f;
pass(13) = normest(g - 1) < 10*epslevel(f);

f = chebfun(@(x) sin(x));
g = (2i).^f;
h = chebfun(@(x) (2i).^sin(x));
pass(14) = normest(g - h) < 10*epslevel(f);

f = chebfun(@(x) [sin(x), cos(x), 1i*exp(x)]);
fq = quasimatrix(f);
g = 1.^f;
pass(15) = min(size(g)) == 3 && normest(g - 1) < 10*epslevel(h);
gq = 1.^fq;
pass(16) = normest(gq - g) < 10*epslevel(h)*vscale(h);

g = (2i).^f;
h = chebfun(@(x) [2i.^sin(x), 2i.^cos(x), 2i.^(1i*exp(x))]);
pass(17) = min(size(g)) == 3 && normest(g - h) < 100*epslevel(h);
gq = (2i).^fq;
pass(18) = normest(gq - g) < 10*epslevel(h)*vscale(h);

%% CHEBFUN .^ CHEBFUN

x = chebfun(@(x) x, [.1, 2]);
g = x.^x;
h = chebfun(@(x) x.^x, [.1, 2]);
pass(19) = normest(g - h, [.1, 2]) < 10*epslevel(h)*vscale(h);

x = chebfun(@(x) [x, exp(1i*x)], [.1, 2]);
g = x.^x;
h = chebfun(@(x) [x.^x, exp(1i*x).^exp(1i*x)], [.1, 2]);
pass(20) = normest(g - h, [.1, 2]) < 10*epslevel(h)*vscale(h);

xq = quasimatrix(x);
gq = xq.^xq;
pass(21) = normest(g - gq, [.1, 2]) < 10*epslevel(h)*vscale(h);

%%%%%%%%%%%%%%%%%%%%%%%%% Integration of SINGFUN: %%%%%%%%%%%%%%%%%%%%%%%%%%
dom = [-2 7];

% Generate a few random points to use as test values:
seedRNG(6178);
domCheck = [dom(1)+0.1 dom(2)-0.1];
x = diff(domCheck) * rand(100, 1) + domCheck(1);

%% SINGFUN .^ 2:
pow = -1.5;
b = 2;
op = @(x) (x-dom(1)).^pow;
opExact = @(x) (x-dom(1)).^(b*pow);

pref.singPrefs.exponents = [pow 0];
f = chebfun(op, dom, pref);
g = f.^b;
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = vals_g - vals_exact;
pass(22) = ( norm(err, inf) < 1e1*epslevel(f).*norm(vals_exact, inf) );

%% SINGFUN .^ constant integer:

pow = -0.5;
b = 3;
op = @(x) sin(100*x).*(x-dom(1)).^pow;
opExact = @(x) sin(100*x).^b.*(x-dom(1)).^(b*pow);

pref.singPrefs.exponents = [pow 0];
pref.enableBreakpointDetection = 1;
f = chebfun(op, dom, pref);
g = f.^b;
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = vals_g - vals_exact;
pass(23) = ( norm(err, inf) < epslevel(f).*norm(vals_exact, inf) );

%% Square root via POWER - A positive piece-wise example with singularities:

pow = -1.5;
op = @(x) (sin(100*x).^2+1).*(x-dom(1)).^pow;
opExact = @(x) sqrt(sin(100*x).^2+1).*(x-dom(1)).^(pow/2);

pref.singPrefs.exponents = [pow 0];
pref.enableBreakpointDetection = 1;
f = chebfun(op, dom, pref);
g = f.^(1/2);
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = vals_g - vals_exact;
pass(24) = ( norm(err, inf) < 1e2*epslevel(f).*norm(vals_exact, inf) );

%% General power - A smooth function - the Runge function and positive power:
% Define a domain:
dom = [-2 7];

% Generate a few random points to use as test values:
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

op = @(x) 1./(1+25*x.^2);
b = 0.6;
opExact = @(x) 1./((1+25*x.^2).^b);

f = chebfun(op, dom, 'splitting', 'on');
g = f.^b;
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = vals_g - vals_exact;
pass(25) = ( norm(err, inf) < epslevel(f).*norm(vals_exact, inf) );

%% General power - A smooth function - the Runge function and negative power:

op = @(x) 1./(1+25*x.^2);
b = -1.6;
opExact = @(x) 1./((1+25*x.^2).^b);

f = chebfun(op, dom, 'splitting', 'on');
g = f.^b;
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = vals_g - vals_exact;
pass(26) = ( norm(err, inf) < 1e2*epslevel(f).*norm(vals_exact, inf) );

%% General power - A smooth function - a real function with varying sign and 
% positive power:

op = @(x) sin(10*x);
b = 0.8;
opExact = @(x) sin(10*x).^b;

f = chebfun(op, dom, 'splitting', 'on');
g = f.^b;
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = norm(vals_g - vals_exact, inf);
tol = 1e1*epslevel(f).*norm(vals_exact, inf);
pass(27) = ( err < tol );

%% General power - A smooth function - a real function with varying sign and 
% negative power:

op = @(x) sin(10*x);
b = -0.65;
opExact = @(x) sin(10*x).^b;

f = chebfun(op, dom, 'splitting', 'on');
g = f.^b;
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = norm(vals_g - vals_exact, inf);
tol = 1e3*epslevel(f).*norm(vals_exact, inf);
pass(28) = ( err < tol );

%% General power - A smooth function - a complex function and positive power:

op = @(x) sin(3*x) + 1i*cos(2*x);
b = 1.7;
opExact = @(x) (sin(3*x) + 1i*cos(2*x)).^b;

f = chebfun(op, dom, 'splitting', 'on');
g = f.^b;
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = vals_g - vals_exact;
pass(29) = ( norm(err, inf) < 1e1*epslevel(f).*norm(vals_exact, inf) );

%% General power - A smooth function - a complex function and negative power:

op = @(x) sin(3*x) + 1i*cos(2*x);
b = -0.77;
opExact = @(x) (sin(3*x) + 1i*cos(2*x)).^b;

f = chebfun(op, dom, 'splitting', 'on');
g = f.^b;
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = vals_g - vals_exact;
pass(30) = ( norm(err, inf) < 1e1*epslevel(f).*norm(vals_exact, inf) );

%% General power - A singular function - a positive function and positive power:
dom = [-2 7];

% Generate a few random points to use as test values:
seedRNG(6178);
domCheck = [dom(1)+0.1 dom(2)-0.1];
x = diff(domCheck) * rand(100, 1) + domCheck(1);

b = 0.77;
op = @(x) ((x-dom(1)).^-1.2);
opExact = @(x) ((x-dom(1)).^(-1.2*b));

f = chebfun(op, dom, 'splitting', 'on', 'exps', [-1.2 0]);
g = f.^b;
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = vals_g - vals_exact;
pass(31) = ( norm(err, inf) < 1e1*epslevel(f).*norm(vals_exact, inf) );

%% General power - A smooth function with varying sign and positive power:

b = 0.77;
op = @(x) sin(20*x).*((x-dom(1)).^-1.2);
opExact = @(x) (sin(20*x).^b).*((x-dom(1)).^(-1.2*b));

f = chebfun(op, dom, 'splitting', 'on', 'exps', [-1.2 0]);
g = f.^b;
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = norm(vals_g - vals_exact, inf);
tol = 1e2*epslevel(f).*norm(vals_exact, inf);
pass(32) = ( err < tol );

%% General power - A smooth function with varying sign and positive power:

b = 0.69;
op = @(x) (sin(20*x)+1i*cos(3*x)).*((x-dom(1)).^-1.2).*((dom(2)-x).^-0.49);
opExact = @(x) ((sin(20*x)+1i*cos(3*x)).^b).*((x-dom(1)).^(-1.2*b)).* ...
    ((dom(2)-x).^(-0.49*b));

f = chebfun(op, dom, 'splitting', 'on', 'exps', [-1.2 -0.49]);
g = f.^b;
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = vals_g - vals_exact;
pass(33) = ( norm(err, inf) < 1e1*epslevel(f).*norm(vals_exact, inf) );

end

function out = normest(f, dom)

% Generate a few random points to use as test values.
seedRNG(6178);
if ( nargin == 1 )
    x = 2 * rand(100, 1) - 1;
else
    x = sum(dom) * rand(10, 1) - dom(1);
end

out = norm(feval(f, x), inf);

end

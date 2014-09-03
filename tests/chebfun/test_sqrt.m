function pass = test_sqrt(pref)

if ( nargin == 0 ) 
    pref = chebfunpref();
end

% Define a domain:
dom = [-2 7];

% Generate a few random points to use as test values:
seedRNG(6178);
x = sort(diff(dom) * rand(100, 1) + dom(1));

%% A positive function (the Runge function):
op = @(x) 1./(1+25*x.^2);
opExact = @(x) 1./sqrt(1+25*x.^2);

f = chebfun(op, dom, pref);
g = sqrt(f);
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = vals_g - vals_exact;
pass(1) = ( norm(err, inf) < 1e2*epslevel(f).*norm(vals_exact, inf) );

%% A piece-wise smooth case:
op = @(x) sin(50*x).^2+1;
opExact = @(x) sqrt(sin(50*x).^2+1);

f = chebfun(op, dom, 'splitting', 'on');
g = sqrt(f);
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = vals_g - vals_exact;
pass(2) = ( norm(err, inf) < 1e2*epslevel(f).*norm(vals_exact, inf) );

%% A complex piece-wise case:

op = @(x) sin(50*x);
opExact = @(x) sqrt(sin(50*x));

pref.splitting = 1;
f = chebfun(op, dom, pref);
g = sqrt(f);
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = norm(vals_g - vals_exact, inf);
tol = 1e1*epslevel(f).*norm(vals_exact, inf);
pass(3) = err < tol;

%% Another complex piece-wise case:

op = @(x) sin(50*x)+1i*cos(30*x);
opExact = @(x) sqrt(sin(50*x)+1i*cos(30*x));

pref.splitting = 1;
f = chebfun(op, dom, pref);
g = sqrt(f);
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = vals_g - vals_exact;
pass(4) = ( norm(err, inf) < 1e2*epslevel(f).*norm(vals_exact, inf) );

%% An array-valued CHEBFUN: 

op = @(x) [sin(x), sin(x)-.5];
opExact = @(x) [sqrt(sin(x)), sqrt(sin(x)-.5)];
f = chebfun(op, dom);
g = sqrt(f);
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = norm(vals_g - vals_exact, inf);
tol = 1e2*get(f,'epslevel')*norm(vals_exact, inf);
pass(5) = err < tol;

%% A positive piece-wise example with singularities:

% Generate a few random points to use as test values:
domCheck = [dom(1)+0.1 dom(2)-0.1];
x = diff(domCheck) * rand(100, 1) + domCheck(1);

pow = -1.5;
op = @(x) (sin(50*x).^2+1).*(x-dom(1)).^pow;
opExact = @(x) sqrt(sin(50*x).^2+1).*(x-dom(1)).^(pow/2);

f = chebfun(op, dom, 'exps', [pow 0], 'splitting', 'on');
g = sqrt(f);
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = vals_g - vals_exact;
pass(6) = ( norm(err, inf) < 1e2*epslevel(f).*norm(vals_exact, inf) );

%%%%%%%%%%%%%%%%%%%%%%% function on unbounded domain: %%%%%%%%%%%%%%%%%%%%%%%%%%

%% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];
domCheck = [-1e2 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

% A function converging to constant values at -Inf and Inf:
op = @(x) x.^2.*exp(-x.^2)+2;
opg = @(x) sqrt(x.^2.*exp(-x.^2)+2);
f = chebfun(op, dom);
g = sqrt(f);
gVals = feval(g, x);
gExact = opg(x);
err = gVals - gExact;
pass(7) = norm(err, inf) < epslevel(g)*vscale(g);

% Blow-up function:
op = @(x) x.^2.*(1-exp(-x.^2))+2;
opg = @(x) sqrt(x.^2.*(1-exp(-x.^2))+2);
f = chebfun(op, dom, 'exps', [2 2]);
g = sqrt(f);
gVals = feval(g, x);
gExact = opg(x);
err = gVals - gExact;
pass(8) = norm(err, inf) < 1e2*epslevel(g)*vscale(g);

%% Functions on [a inf]:

% Set the domain:
dom = [1 Inf];
domCheck = [1 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

% Integer power of a function converging to a constant at Inf:
op = @(x) x.*exp(-x)+3;
opg = @(x) sqrt(x.*exp(-x)+3);
f = chebfun(op, dom);
g = sqrt(f);
gVals = feval(g, x);
gExact = opg(x);
err = gVals - gExact;
pass(9) = norm(err, inf) < epslevel(g)*vscale(g);

% Oscillatory function with varying sign and integer power:
op = @(x) 0.1+sin(10*x)./exp(x);
opg = @(x) sqrt(0.1+sin(10*x)./exp(x));
f = chebfun(op, dom, 'splitting', 'on');
g = sqrt(f);
gVals = feval(g, x);
gExact = opg(x);
err = gVals - gExact;
pass(10) = norm(err, inf) < epslevel(g)*vscale(g);

% Blow-up function and negative integer power:
op = @(x) x.*(5+exp(-x.^3));
opg = @(x) sqrt(x.*(5+exp(-x.^3)));
f = chebfun(op, dom, 'exps', [0 1]);
g = sqrt(f);
gVals = feval(g, x);
gExact = opg(x);
err = gVals - gExact;
pass(11) = norm(err, inf) < epslevel(g)*vscale(g);

end

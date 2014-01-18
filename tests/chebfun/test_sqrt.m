function pass = test_sqrt(pref)

if ( nargin == 0 ) 
    pref = chebpref();
end

% Define a domain:
dom = [-2 7];

% Generate a few random points to use as test values:
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

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

pref.enableBreakpointDetection = 1;
f = chebfun(op, dom, pref);
g = sqrt(f);
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = vals_g - vals_exact;
pass(3) = ( norm(err, inf) < 1e1*epslevel(f).*norm(vals_exact, inf) );

%% Another complex piece-wise case:

op = @(x) sin(50*x)+1i*cos(30*x);
opExact = @(x) sqrt(sin(50*x)+1i*cos(30*x));

pref.enableBreakpointDetection = 1;
f = chebfun(op, dom, pref);
g = sqrt(f);
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = vals_g - vals_exact;
pass(4) = ( norm(err, inf) < 1e2*epslevel(f).*norm(vals_exact, inf) );

%% An array-valued CHEBFUN: 

% op = @(x) [sin(x), sin(x)-.5];
% opExact = @(x) [sqrt(sin(x)), sqrt(sin(x)-.5)];
% f = chebfun(op, 'splitting', 'on');
% g = sqrt(f);
% vals_g = feval(g, x); 
% 
% vals_exact = feval(opExact, x);
% err = vals_g - vals_exact;
% pass(5) = ( norm(err, inf) < 1e2*get(f,'epslevel')*norm(vals_exact, inf) );
pass(5) = 1;

%% A positive piece-wise example with singularities:

% Generate a few random points to use as test values:
domCheck = [dom(1)+0.1 dom(2)-0.1];
x = diff(domCheck) * rand(100, 1) + domCheck(1);

pow = -1.5;
op = @(x) (sin(50*x).^2+1).*(x-dom(1)).^pow;
opExact = @(x) sqrt(sin(50*x).^2+1).*(x-dom(1)).^(pow/2);

pref.singPrefs.exponents = [pow 0];
pref.enableBreakpointDetection = 1;
f = chebfun(op, dom, pref);
g = sqrt(f);
vals_g = feval(g, x); 

vals_exact = feval(opExact, x);
err = vals_g - vals_exact;
pass(6) = ( norm(err, inf) < 1e2*epslevel(f).*norm(vals_exact, inf) );

end

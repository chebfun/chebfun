% Test file for @classicfun/rdivide.m

function pass = test_rdivide(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

singPref = pref;
singPref.blowup = true;

% Set a domain for BNDFUN.
data.domain = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(data.domain) * rand(1000, 1) + data.domain(1);

% Random numbers to use as arbitrary constants.
alpha = -0.194758928283640 + 0.075474485412665i;
beta = -0.526634844879922 - 0.685484380523668i;

%%
% Check division by single scalars.
f_op = @(x) sin(x);
f = bndfun(f_op, data, pref);
pass(1) = test_div_function_by_scalar(f, f_op, alpha, x);
    
g = f ./ 0;
pass(2) = isnan(g);
    
f_op = @(x) [sin(x) cos(x)];
f = bndfun(f_op, data, pref);
pass(3) = test_div_function_by_scalar(f, f_op, alpha, x);
    
g = f ./ 0;
pass(4) = isnan(g);
    
%%
% Check division by a row matrix of scalars in the case of an array-valued
% BNDFUN object.
g = f ./ [alpha beta];
g_exact = @(x) [sin(x)./alpha cos(x)./beta];
pass(5) = norm(feval(g, x) - g_exact(x), inf) < ...
    10*max(get(g, 'vscale').*get(g, 'epslevel'));
    
g = f ./ [alpha 0];
isn = isnan(feval(g, x));
pass(6) = isnan(g) && ~any(any(isn(:,1))) ...
    && all(isn(:,2));
    
%%
% Check division of a scalar by a BNDFUN object.
f_op = @(x) exp(x);
f = bndfun(@(x) exp(x), data, pref);
pass(7) = test_div_scalar_by_function(alpha, f, f_op, x);
    
%%
% Check division of two BNDFUN objects.
g_op = @(x) exp(x);
g = bndfun(g_op, data, pref);
    
f_op = @(x) exp(x) - 1;
f = bndfun(f_op, data, pref);
pass(8) = test_div_function_by_function(f, f_op, g, g_op, x);
    
f_op = @(x) 1./(1 + x.^2);
f = bndfun(f_op, data, pref);
pass(9) = test_div_function_by_function(f, f_op, g, g_op, x);
    
f_op = @(x) cos(1e4*x);
f = bndfun(f_op, data, pref);
pass(10) = test_div_function_by_function(f, f_op, g, g_op, x);
    
f_op = @(t) sinh(t*exp(2*pi*1i/6));
f = bndfun(f_op, data, pref);
pass(11) = test_div_function_by_function(f, f_op, g, g_op, x);
    
%%
% Check that direct construction and RDIVIDE give comparable results.
f = bndfun(@(x) sin(x), data, pref);
h1 = f ./ alpha;
h2 = bndfun(@(x) sin(x) ./ alpha, data, pref);
pass(12) = norm(feval(h1, x) - feval(h2, x), inf) < ...
    10*get(h2, 'vscale')*get(h2, 'epslevel');
    
g = bndfun(@(x) exp(x), data, pref);
h1 = f ./ g;
h2 = bndfun(@(x) sin(x) ./ exp(x), data, pref);
pass(13) = norm(feval(h1, x) - feval(h2, x), inf) < ...
    5e3*get(h2, 'vscale')*get(h2, 'epslevel');

%% 
% Test on singular BNDFUN.
pow1 = -0.5;
pow2 = -0.3;
op1 = @(x) (x - data.domain(2)).^pow1.*sin(x);
op2 = @(x) (x - data.domain(2)).^pow2.*(cos(x).^2+1);
singData = data;
singData.exponents = [0 pow1];
f = bndfun(op1, singData, singPref);
singData.exponents = [0 pow2];
g = bndfun(op2, singData, singPref);
h = f./g;
vals_h = feval(h, x);
pow = pow1-pow2;
op = @(x)  (x - data.domain(2)).^pow.*(sin(x)./(cos(x).^2+1));
h_exact = op(x);
pass(14) = ( norm(vals_h-h_exact, inf) < 1e2* ...
    max(get(f, 'epslevel'), get(g, 'epslevel'))*norm(h_exact, inf) );

%% Tests for UNBNDFUN:

% Functions on [2 inf]:

% Set the domain:
data.domain = [2 Inf];
domCheck = [2 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

opf = @(x) exp(-x.^2);
opg = @(x) x.^2;
oph = @(x) exp(-x.^2).*x.^-2;
f = unbndfun(opf, data);
g = unbndfun(opg, data);
h = f./g;
hVals = feval(h, x);
hExact = oph(x);
err = hVals - hExact;
pass(15) = norm(err, inf) < 1e1*get(f,'epslevel')*get(f,'vscale');

end

%%
% Test the division of a BNDFUN F, specified by F_OP, by a scalar ALPHA using
% a grid of points X in [-1  1] for testing samples.
function result = test_div_function_by_scalar(f, f_op, alpha, x)
    g = f ./ alpha;
    g_exact = @(x) f_op(x) ./ alpha;
    err = norm(feval(g, x) - g_exact(x), inf);
    tol = 10*max(get(g, 'vscale').*get(g, 'epslevel'));
    result = err < tol;
end

%%
% Test the division of a scalar ALPHA by a BNDFUN, specified by F_OP, using
% a grid of points X in [-1  1] for testing samples.
function result = test_div_scalar_by_function(alpha, f, f_op, x)
    g = alpha ./ f;
    g_exact = @(x) alpha ./ f_op(x);
    err = norm(feval(g, x) - g_exact(x), inf);
    tol = 300*max(get(g, 'vscale').*get(g, 'epslevel'));
    result = err < tol;
end

%%
% Test the division of two BNDFUN objects F and G, specified by F_OP and
% G_OP, using a grid of points X in [-1  1] for testing samples.
function result = test_div_function_by_function(f, f_op, g, g_op, x)
    h = f ./ g;
    h_exact = @(x) f_op(x) ./ g_op(x);
    err = norm(feval(h, x) - h_exact(x), inf);
    tol = 300*max(get(h, 'vscale').*get(h, 'epslevel'));
    result = err < tol;
end

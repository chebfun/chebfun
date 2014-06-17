% Test file for bndfun/feval.m

function pass = test_feval(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Set the domain
dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(7681);
x = diff(dom) * rand(1000, 1) + dom(1);

%%
% Spot-check values for a couple of functions.  We can only expect
% accuracy on the order of the truncation level, so we use this as our
% criterion.

f = bndfun(@(x) exp(x) - 1, struct('domain', dom), pref);
f_exact = @(x) exp(x) - 1;
pass(1) = (norm(feval(f, x) - f_exact(x), inf) < ...
    2*get(f, 'vscale')*get(f, 'epslevel'));

f = bndfun(@(x) 1./(1 + x.^2), struct('domain', dom), pref);
f_exact = @(x) 1./(1 + x.^2);
pass(2) = (norm(feval(f, x) - f_exact(x), inf) < ...
    get(f, 'vscale')*get(f, 'epslevel'));

f = bndfun(@(x) cos(1e4*x), struct('domain', dom), pref);
f_exact = @(x) cos(1e4*x);
pass(3) = (norm(feval(f, x) - f_exact(x), inf) < 1e4*get(f, 'epslevel'));

z = exp(2*pi*1i/6);
f = bndfun(@(t) sinh(t*z), struct('domain', dom), pref);
f_exact = @(t) sinh(t*z);
pass(4) = (norm(feval(f, x) - f_exact(x), inf) < ...
    10*get(f, 'vscale')*get(f, 'epslevel'));

%%
% Check row vector and matrix input.

err = feval(f, x.') - f_exact(x.');
pass(5) = (all(size(err) == [1 1000])) && (norm(err(:), inf) < ...
    10*get(f, 'vscale')*get(f, 'epslevel'));

x_mtx = reshape(x, [100 10]);
err = feval(f, x_mtx) - f_exact(x_mtx);
pass(6) = (all(size(err) == [100 10])) && (norm(err(:), inf) < ...
    10*get(f, 'vscale')*get(f, 'epslevel'));

x_3mtx = reshape(x, [10 10 10]);
err = feval(f, x_3mtx) - f_exact(x_3mtx);
pass(7) = (all(size(err) == [10 10 10])) && (norm(err(:), inf) < ...
    10*get(f, 'vscale')*get(f, 'epslevel'));

%%
% Check operation for array-valued bndfun objects.

f = bndfun(@(x) [sin(x) x.^2 exp(1i*x)], struct('domain', dom), pref);
f_exact = @(x) [sin(x) x.^2 exp(1i*x)];
err = feval(f, x) - f_exact(x);
pass(8) = all(max(abs(err)) < 10*get(f, 'vscale').*get(f, 'epslevel'));

%%
% Test for evaluating array-valued bndfun objects at matrix arguments if
% the operation makes sense.

f = bndfun(@(x) [sin(pi*x) cos(pi*x)], struct('domain', dom), pref);
x2 = [-1 0 5 ; -1.75 .5 4.75];
fx = feval(f, x2);
f_exact = [0 0 0 -1 1 -1
    [1 sqrt(2) 1 1 0 -1]/sqrt(2)];
pass(9) = all(all(abs(fx - f_exact) < ...
    10*max(get(f, 'vscale').*get(f, 'epslevel'))));

%% Test on singular function:

pow = -0.5;
op = @(x) (x - dom(1)).^pow.*sin(x);
pref.blowup = true;
data.domain = dom;
data.exponents = [pow 0];
f = bndfun(op, data, pref);
fval = feval(f, x);
vals_exact = feval(op, x);
err = fval - vals_exact;
pass(10) = ( norm(err, inf) < 1e1*get(f,'epslevel')*norm(vals_exact, inf) );

end

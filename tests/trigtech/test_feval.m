% Test file for trigtech/feval.m

function pass = test_feval(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

% Generate a few random points to use as test values.
seedRNG(7681);
x = 2 * rand(1000, 1) - 1;

testclass = trigtech();

%%
% Spot-check values for a couple of functions.  We can only expect 
% accuracy on the order of the truncation level, so we use this as our 
% criterion.

f = testclass.make(@(x) sin(pi*x), [], pref);
f_exact = @(x) sin(pi*x);
pass(1) = (norm(feval(f, x) - f_exact(x), inf) < ...
    10*f.vscale.*f.epslevel);

f = testclass.make(@(x) exp(cos(pi*x)) - 1, [], pref);
f_exact = @(x) exp(cos(pi*x)) - 1;
pass(2) = (norm(feval(f, x) - f_exact(x), inf) < ...
    10*f.vscale.*f.epslevel);

f = testclass.make(@(x) cos(100*sin(pi*x)), [], pref);
f_exact = @(x) cos(100*sin(pi*x));
pass(3) = (norm(feval(f, x) - f_exact(x), inf) < ...
    10*f.vscale.*f.epslevel);

f = testclass.make(@(x) exp(1i*pi*x), [], pref);
f_exact = @(x) exp(1i*pi*x);
pass(4) = (norm(feval(f, x) - f_exact(x), inf) < ...
    10*f.vscale.*f.epslevel);
    
%%
% Check row vector and matrix input.
    
err = feval(f, x.') - f_exact(x.');
pass(5) = (all(size(err) == [1 1000])) && (norm(err(:), inf) < ...
    10*f.vscale.*f.epslevel);

x_mtx = reshape(x, [100 10]);
err = feval(f, x_mtx) - f_exact(x_mtx);
pass(6) = (all(size(err) == [100 10])) && (norm(err(:), inf) < ...
    10*f.vscale.*f.epslevel);

x_3mtx = reshape(x, [10 10 10]);
err = feval(f, x_3mtx) - f_exact(x_3mtx);
pass(7) = (all(size(err) == [10 10 10])) && (norm(err(:), inf) < ...
    10*f.vscale.*f.epslevel);

%%
% Check operation for array-valued trigtech objects.

f = testclass.make(@(x) [(2+sin(pi*x)).*exp(1i*pi*x), -(2+sin(pi*x)).*exp(1i*pi*x), 2+sin(pi*x)], [], pref);
f_exact = @(x) [(2+sin(pi*x)).*exp(1i*pi*x), -(2+sin(pi*x)).*exp(1i*pi*x), 2+sin(pi*x)];
err = feval(f, x) - f_exact(x);
pass(8) = all(max(abs(err)) < 10*max(f.vscale.*f.epslevel));

%%
% Test for evaluating array-valued trigtech objects at matrix arguments if
% the operation makes sense.

f = testclass.make(@(x) [sin(pi*x) cos(pi*x) exp(1i*pi*x)], [], pref);
x2 = [-1 0 1 ; .25 .5 .75];
fx = feval(f, x2);
f_exact = [0 0 0 -1 1 -1 exp(-1i*pi) 1 exp(1i*pi)
    [1 sqrt(2) 1 1 0 -1]/sqrt(2) exp(1i*pi.*[.25 .5 .75])];
pass(9) = all(all(abs(fx - f_exact) < 10*max(f.vscale.*f.epslevel)));

end

% Test file for singfun/rdivide.m

function pass = test_rdivide(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebpref();
end

% Generate a few random points to use as test values.
seedRNG(666);
x = 2 * rand(100, 1) - 1;
x = sort(x);

%%
% Check operation in the case of empty arguments.
f = singfun();
g = singfun(@(x) 1./(1+x), [-1, 0]);
pass(1) = (isempty(f./g) && isempty(g./f) && isempty(f./f));

%%
% SMOOTHFUN ./ SINGFUN
f = smoothfun.constructor(@(x) 2 + sin(x));
g = singfun(@(x) cos(x));
pass(2) = isa(f./g, 'smoothfun') && isa(g./f, 'smoothfun');
%%
% Division of smooth SINGFUNs should not return a SINGFUN
f = singfun(@(x) sin(x));
g = singfun(@(x) cos(x));
pass(3) = ~isa(f./g, 'singfun');
%%
% Check division with a double.
fh = @(x) 1./((1+x).*(1-x));
f = singfun(fh, [-1, -1]);
% A random double:
g = rand();
pass(4) = test_division_by_scalar(f, fh, g, x);

%%
% Check reciprocal of a singfun.
fh = @(x) 0*x + 1;
f = singfun(fh, [], {'none', 'none'}, [], [], pref);
gh = @(x) cos(x);
g = singfun(gh, [], {'none', 'none'}, [], [], pref);
pass(5) = test_divide_function_by_function(f, fh, g, gh, x);

%% 
% Check division of two singfun objects.
fh = @(x) cos(x);
f = singfun(fh, [], {'none', 'none'}, [], [], pref);
pass(6) = test_divide_function_by_function(f, fh, f, fh, x);

fh = @(x) sin(x);
f = singfun(fh, [], [], [], [], pref);

gh = @(x) (1+x).*(1-x);  %
g = singfun(gh, [], [], [], [], pref);
% Remove points really close to the end points.
pass(7) = test_divide_function_by_function(f, fh, g, gh, x(5:end-4));

fh = @(x) sin(1e2*x);
f = singfun(fh, [], [], [], [], pref);
% Remove points really close to the end points.
pass(8) = test_divide_function_by_function(f, fh, g, gh, x(5:end-4));

%%
% Check that direct construction and RDIVIDE give comparable results.

f = singfun(@(x) sin(x), [], [], [], [], pref);
g = singfun(@(x) (1+x).*(1-x), [], [], [], [], pref);
h1 = f./g;
h2 = singfun(@(x) sin(x)./((1+x).*(1-x)), [], [], [], [], pref);
x = x(5:end-4); % Remove some points close to the end points.
tol = 1e2*max(get(h1, 'vscale')*get(h1, 'epslevel'), ...
    get(h2, 'vscale')*get(h2, 'epslevel'));
pass(9) = norm(feval(h1, x) - feval(h2, x), inf) < tol;


end

% Test the division of a SINGFUN F, specified by Fh, by a scalar C using
% a grid of points X for testing samples.
function result = test_division_by_scalar(f, fh, c, x)
    g = f./c;
    g_exact = @(x) fh(x)./c;
    result = norm(feval(g, x) - g_exact(x), inf) <= ...
        1e2*get(g, 'vscale')*get(g, 'epslevel');
end

% Test the division of two SINGFUN objects F and G, specified by FH and
% GH, using a grid of points X for testing samples.
function result = test_divide_function_by_function(f, fh, g, gh, x)
    h = f./g;
    h_exact = @(x) fh(x)./gh(x);
    result = norm(feval(h, x) - h_exact(x), inf) <= ...
        1e2*get(h, 'vscale')*get(h, 'epslevel');
end

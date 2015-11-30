% Test file for singfun/rdivide.m

function pass = test_rdivide(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(666);
x = 2 * rand(100, 1) - 1;
x = sort(x);

%%
% Check operation in the case of empty arguments.
f = singfun();
data.exponents = [-1 0];
g = singfun(@(x) 1./(1+x), data);
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
data.exponents = [-1 -1];
f = singfun(fh, data);
% A random double:
g = rand();
pass(4) = test_division_by_scalar(f, fh, g, x);

%%
% Check reciprocal of a singfun.
fh = @(x) 0*x + 1;
data.exponents = [];
data.singType = {'none', 'none'};
f = singfun(fh, data, pref);
gh = @(x) cos(x);
data.exponents = [];
data.singType = {'none', 'none'};
g = singfun(gh, data, pref);
pass(5) = test_divide_function_by_function(f, fh, g, gh, x);

%% 
% Check division of two singfun objects.
fh = @(x) cos(x);
data.exponents = [];
data.singType = {'none', 'none'};
f = singfun(fh, data, pref);
pass(6) = test_divide_function_by_function(f, fh, f, fh, x);

fh = @(x) sin(x);
f = singfun(fh, [], pref);

gh = @(x) (1+x).*(1-x);  %
g = singfun(gh, [], pref);
% Remove points really close to the end points.
pass(7) = test_divide_function_by_function(f, fh, g, gh, x(5:end-4));

fh = @(x) sin(1e2*x);
f = singfun(fh, [], pref);
% Remove points really close to the end points.
pass(8) = test_divide_function_by_function(f, fh, g, gh, x(5:end-4));
    
%%
% Check that direct construction and RDIVIDE give comparable results.

f = singfun(@(x) sin(x), [], pref);
g = singfun(@(x) (1+x).*(1-x), [], pref);
h1 = f./g;
h2 = singfun(@(x) sin(x)./((1+x).*(1-x)), [], pref);
x = x(5:end-4); % Remove some points close to the end points.
tol = 1e2*max(get(h1, 'vscale')*eps, ...
    get(h2, 'vscale')*eps);
pass(9) = norm(feval(h1, x) - feval(h2, x), inf) < tol;

%%
% Check exponents of reciprocals of singfuns.
a = randi(20); b = randi(20);

% flip the exponents from positive to negative
gh = @(x) ((1+x).^a).*((1-x).^b);
data.exponents = [a b];
data.singType = {'pole', 'pole'};
g = singfun(gh, data, pref);
h = 1./g;
pass(10) = all( h.exponents == [-a -b] );

% simplify / smooth out exponents when greater than or equal to 1
gh = @(x) ((1+x).^-a).*((1-x).^-b);
data.exponents = [];
data.singType = {'pole', 'pole'};
g = singfun(gh, data, pref);
h = 1./g;
if isa(h,'singfun')
    pass(11) = all( get(h, 'exponents') < 1);
else
    pass(11) = true;
end

%% Check division as differentiation
dom = [1, Inf]; 
r = chebfun(@(x) x, dom); 
f = 1./r;
err = diff(f) + f./r;
pass(12) = ( norm(err, Inf) < 1e3*get(r, 'vscale')*eps );
    

%% Check division as negative powers
f = 1./r - r.^-1;               % rdivide == power
g = (1./r)./r - r.^-2;          % rdivides == times
h = 1./(r.^2) - (1./r).^2;      % associativity
pass(13) = ( norm(f, Inf) < get(r, 'vscale')*eps );
pass(14) = ( norm(g, Inf) < get(r, 'vscale')*eps );
pass(15) = ( norm(h, Inf) < 1e2*get(r, 'vscale')*eps );
    
%% Check simplification of (1+x) / sqrt(1+x) has positive exponents.
f = singfun(@(x) 1+x);
g = singfun(@(x) sqrt(1+x), struct('exponents', [.5 0]));
h = f./g;
pass(16) = all(h.exponents == [.5 0]);

end

% Test the division of a SINGFUN F, specified by Fh, by a scalar C using
% a grid of points X for testing samples.
function result = test_division_by_scalar(f, fh, c, x)
    g = f./c;
    g_exact = @(x) fh(x)./c;
    result = norm(feval(g, x) - g_exact(x), inf) <= ...
        2e3*get(g, 'vscale')*eps;
end

% Test the division of two SINGFUN objects F and G, specified by FH and
% GH, using a grid of points X for testing samples.
function result = test_divide_function_by_function(f, fh, g, gh, x)
    h = f./g;
    h_exact = @(x) fh(x)./gh(x);
    result = norm(feval(h, x) - h_exact(x), inf) <= ...
        1e3*get(h, 'vscale')*eps;
end

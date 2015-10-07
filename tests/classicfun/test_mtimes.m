% Test file for @classicfun/mtimes.m

function pass = test_mtimes(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Set a domain for BNDFUN.
data.domain = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(data.domain) * rand(1000, 1) + data.domain(1);

% A random number to use as an arbitrary scalar multiplier.
alpha = randn() + 1i*randn();

%%
% Check operation in the face of empty BNDFUN.
f = bndfun(@(x) sin(x), data, pref);
g = bndfun();
pass(1) = isempty(f*[]) && isempty([]*f) && isempty(2*g) && isempty(g*2);
    
%%
% Check operation for scalar BNDFUN objects.    
f = bndfun(@(x) sin(x), data, pref);
g1 = alpha*f;
g2 = f*alpha;
pass(2) = isequal(g1, g2);
g_exact = @(x) alpha*sin(x);
pass(3) = norm(feval(g1, x) - g_exact(x), inf) < ...
    10*get(g1, 'vscale')*eps;

g = 0*f;
pass(4) = all(feval(g, x) == 0);
    
%%
% Check operation for array-valued BNDFUN objects.
f = bndfun(@(x) [sin(x) cos(x) exp(x)], data, pref);
g1 = alpha*f;
g2 = f*alpha;
pass(5) = isequal(g1, g2);
g_exact = @(x) alpha*[sin(x) cos(x) exp(x)];
err = abs(feval(g1, x) - g_exact(x));
pass(6) = max(err(:)) < 10*max(get(g1, 'vscale')*eps);
    
g = 0*f;
pass(7) = all(all(feval(g, x) == zeros(numel(x), 3)));
    
A = randn(3, 3);
g = f*A;
g_exact = @(x) [sin(x) cos(x) exp(x)]*A;
err = abs(feval(g, x) - g_exact(x));
pass(8) = max(err(:)) < 10*max(get(g, 'vscale')*eps);
    
    
%%
% Verify error handling and corner cases.

% Multiply non-scalar double and fun.
try
    f = bndfun(@(x) exp(x), data);
    disp([1 2 3]*f)
catch ME
    pass(9) = strcmp(ME.identifier, 'CHEBFUN:CLASSICFUN:mtimes:size') ...
        && strcmp(ME.message, 'Inner matrix dimensions must agree.');
end
    
% Multiply fun and non-scalar double with mismatching dimensions.
try
    f = bndfun(@(x) [sin(x) cos(x)], data);
    disp(f*[1 ; 2 ; 3]);
    pass(10) = false;
catch ME
    pass(10) = strcmp(ME.identifier, 'CHEBFUN:CHEBTECH:mtimes:size2') ...
        && strcmp(ME.message, 'Inner matrix dimensions must agree.');
end
    
% Using * for multiplication of two fun objects.
try
    g = bndfun(@(x) x, data);
    disp(f*g);
    pass(11) = false;
catch ME
    pass(11) = strcmp(ME.message, 'Use .* to multiply CLASSICFUN objects.');
end
    
% Using * to multiply a fun and something else.
try
    disp(f*uint8(128));
    pass(12) = false;
catch ME
    pass(12) = strcmp(ME.message, ...
        'mtimes does not know how to multiply a CLASSICFUN and a uint8.');
end

%% Tests for UNBNDFUN:

% Functions on [-inf b]:

% Set the domain:
data.domain = [-Inf -3*pi];
domCheck = [-1e6 -3*pi];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

% Array-valued function:
op = @(x) [exp(x) x.*exp(x) (1-exp(x))./x];
f = unbndfun(op, data);
g = f*A;
gVals = feval(g, x);

op = @(x) [exp(x) x.*exp(x) (1-exp(x))./x]*A;
gExact = op(x);
err = gVals - gExact;
pass(13) = norm(err, inf) < 1e2*max(eps*get(g,'vscale'));


end

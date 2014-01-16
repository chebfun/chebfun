% Test file for fun/mtimes.m

function pass = test_mtimes(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebpref();
end

% Set a domain for BNDFUN.
dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(1000, 1) + dom(1);

% A random number to use as an arbitrary scalar multiplier.
alpha = randn() + 1i*randn();

%%
% Check operation in the face of empty BNDFUN.
f = bndfun(@(x) sin(x), dom, [], [], pref);
g = bndfun();
pass(1) = isempty(f*[]) && isempty([]*f) && isempty(2*g) && isempty(g*2);
    
%%
% Check operation for scalar BNDFUN objects.    
f = bndfun(@(x) sin(x), dom, [], [], pref);
g1 = alpha*f;
g2 = f*alpha;
pass(2) = isequal(g1, g2);
g_exact = @(x) alpha*sin(x);
pass(3) = norm(feval(g1, x) - g_exact(x), inf) < ...
    10*get(g1, 'vscale')*get(g1, 'epslevel');

g = 0*f;
pass(4) = all(feval(g, x) == 0);
    
%%
% Check operation for array-valued BNDFUN objects.
f = bndfun(@(x) [sin(x) cos(x) exp(x)], dom, [], [], pref);
g1 = alpha*f;
g2 = f*alpha;
pass(5) = isequal(g1, g2);
g_exact = @(x) alpha*[sin(x) cos(x) exp(x)];
err = abs(feval(g1, x) - g_exact(x));
pass(6) = max(err(:)) < 10*max(get(g1, 'vscale'))*get(g1, 'epslevel');
    
g = 0*f;
pass(7) = all(all(feval(g, x) == zeros(numel(x), 3)));
    
A = randn(3, 3);
g = f*A;
g_exact = @(x) [sin(x) cos(x) exp(x)]*A;
err = abs(feval(g, x) - g_exact(x));
pass(8) = max(err(:)) < 2*max(get(g, 'vscale'))*get(g, 'epslevel');
    
%%
% Verify error handling and corner cases.
    
% Multiply non-scalar double and fun.
try
    f = bndfun(@(x) exp(x), dom);
    disp([1 2 3]*f)
catch ME
    pass(9) = strcmp(ME.identifier, 'CHEBFUN:FUN:mtimes:size') ...
    && strcmp(ME.message, 'Inner matrix dimensions must agree.');
end
    
% Multiply fun and non-scalar double with mismatching dimensions.
try
    f = bndfun(@(x) [sin(x) cos(x)], dom);        
    disp(f*[1 ; 2 ; 3]);
    pass(10) = false;
catch ME
    pass(10) = strcmp(ME.identifier, 'CHEBFUN:CHEBTECH:mtimes:size2') ...
        && strcmp(ME.message, 'Inner matrix dimensions must agree.');
end
    
% Using * for multiplication of two fun objects.
try
    g = bndfun(@(x) x, dom);
    disp(f*g);
    pass(11) = false;
catch ME
    pass(11) = strcmp(ME.message, 'Use .* to multiply FUN objects.');
end
    
% Using * to multiply a fun and something else.
try
    disp(f*uint8(128));
    pass(12) = false;
catch ME
    pass(12) = strcmp(ME.message, ...
        'mtimes does not know how to multiply a FUN and a uint8.');
end

%% 
% [TODO]: Run a few tests for UNBNDFUN.
end
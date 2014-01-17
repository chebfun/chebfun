% Test file for fun/minus.m

function pass = test_minus(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebpref();
end

% Set a domain for BNDFUN.
dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

% A random number to use as an arbitrary additive constant.
alpha = randn() + 1i*randn();
   
%%
% Check operation in the face of empty arguments.
f = bndfun();
g = bndfun(@(x) x, dom, [], [], pref);
pass(1) = (isempty(f - f) && isempty(f - g) && isempty(g - f));
    
%% 
% Check subtraction with scalars.
f_op = @(x) sin(x);
f = bndfun(f_op, dom, [], [], pref);
pass(2:3) = test_sub_function_and_scalar(f, f_op, alpha, x);
    
%%
% Check subtraction of two BNDFUN objects.
f_op = @(x) zeros(size(x));
f = bndfun(f_op, dom, [], [], pref);
pass(4:5) = test_sub_function_and_function(f, f_op, f, f_op, x);
    
f_op = @(x) exp(x) - 1;
f = bndfun(f_op, dom, [], [], pref);
    
g_op = @(x) 1./(1 + x.^2);
g = bndfun(g_op, dom, [], [], pref);
pass(6:7) = test_sub_function_and_function(f, f_op, g, g_op, x);
    
g_op = @(x) cos(1e4*x);
g = bndfun(g_op, dom, [], [], pref);
pass(8:9) = test_sub_function_and_function(f, f_op, g, g_op, x);
    
g_op = @(t) sinh(t*exp(2*pi*1i/6));
g = bndfun(g_op, dom, [], [], pref);
pass(10:11) = test_sub_function_and_function(f, f_op, g, g_op, x);
    
%% 
% Check operation for array-valued BNDFUN objects.
f_op = @(x) [zeros(size(x)) zeros(size(x)) zeros(size(x))];
f = bndfun(f_op, dom, [], [], pref);
pass(12:13) = test_sub_function_and_function(f, f_op, f, f_op, x);
    
f_op = @(x) [sin(x) cos(x) exp(x)];
f = bndfun(f_op, dom, [], [], pref);
pass(14:15) = test_sub_function_and_scalar(f, f_op, alpha, x);
    
g_op = @(x) [cosh(x) airy(1i*x) sinh(x)];
g = bndfun(g_op, dom, [], [], pref);
pass(16:17) = test_sub_function_and_function(f, f_op, g, g_op, x);
    
%% 
% This should fail with a dimension mismatch error.
g_op = @(x) sin(x);
g = bndfun(g_op, dom, [], [], pref);
try
    h = f - g; %#ok<NASGU>
    pass(18) = false;
catch ME
    pass(18) = strcmp(ME.message, 'Matrix dimensions must agree.');
end
    
%% 
% Check that direct construction and MINUS give comparable results.
f = bndfun(@(x) x, dom, [], [], pref);
g = bndfun(@(x) cos(x) - 1, dom, [], [], pref);
h1 = f - g;
h1_vals = feval(h1, x);
h2 = bndfun(@(x) x - (cos(x) - 1), dom, [], [], pref);
h2_vals = feval(h2, x);
pass(19) = ( norm(h1_vals - h2_vals, inf) < ...
    1e1*get(h1, 'epslevel').*get(h1, 'vscale') );

%% 
% Check that subtracting a BNDFUN and an unhappy BNDFUN gives an
% unhappy result.  
f = bndfun(@(x) cos(x+1), dom);    % Happy
g = bndfun(@(x) sqrt(x+1), dom);   % Unhappy
h = f - g;  % Subtract unhappy from happy.
pass(20) = (~get(g, 'ishappy')) && (~get(h, 'ishappy'));
h = g - f;  % Subtract happy from unhappy.
pass(21) = (~get(g, 'ishappy')) && (~get(h, 'ishappy'));

%% 
% [TODO]: Run a few tests for UNBNDFUN.
end

%% 
% Test the subtraction of a BNDFUN F, specified by F_OP, to and from a scalar
% ALPHA using a grid of points X in [a  b] for testing samples.
function result = test_sub_function_and_scalar(f, f_op, alpha, x)
g1 = f - alpha;
g2 = alpha - f;
result(1) = isequal(g1, -g2);
g_exact = @(x) f_op(x) - alpha;
result(2) = norm(feval(g1, x) - g_exact(x), inf) < ...
    10*max(get(f,'vscale'))*get(g1, 'epslevel');
end

%% 
% Test the subraction of two BNDFUN objects F and G, specified by F_OP and
% G_OP, using a grid of points X in [-1  1] for testing samples.
function result = test_sub_function_and_function(f, f_op, g, g_op, x)
h1 = f - g;
h2 = g - f;
result(1) = isequal(h1, -h2);
h_exact = @(x) f_op(x) - g_op(x);
norm(feval(h1, x) - h_exact(x), inf);
result(2) = norm(feval(h1, x) - h_exact(x), inf) <= 10*max(get(h1,'vscale'))*get(h1,'epslevel');
end
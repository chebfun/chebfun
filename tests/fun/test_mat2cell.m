% Test file for fun/mat2cell.m

function pass = test_mat2cell(pref)

% Get preferences.
if ( nargin < 2 )
    pref = chebpref();
end

% Set a domain for BNDFUN.
dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(1000, 1) + dom(1);

%% 
% Run a few tests for BNDFUN.

f = bndfun(@(x) [sin(x) cos(x) exp(x) x], dom, [], [], pref);
g = bndfun(@(x) sin(x), dom, [], [], pref);
h = bndfun(@(x) [cos(x) exp(x)], dom, [], [], pref);
l = bndfun(@(x) x, dom, [], [], pref);

g_vals = feval(g, x);
h_vals = feval(h, x);
l_vals = feval(l, x);

% Test full arguments.
F = mat2cell(f, 1, [1 2 1]);
F1_vals = feval(F{1}, x);
F2_vals = feval(F{2}, x);
F3_vals = feval(F{3}, x);

pass(1) = ~isempty(F{1}) && ...
    norm(F1_vals - g_vals, inf) < 1e1*get(g, 'epslevel')*get(g, 'vscale');
pass(2) = ~isempty(F{2}) && ...
    norm(F2_vals - h_vals, inf) < max(get(h, 'epslevel')*get(h, 'vscale'));
pass(3) = ~isempty(F{3}) && ...
    norm(F3_vals - l_vals) < 1e1*get(l, 'epslevel')*get(l, 'vscale');
    
% Test two arguments.
F = mat2cell(f, [1 2 1]);
F1_vals = feval(F{1}, x);
F2_vals = feval(F{2}, x);
F3_vals = feval(F{3}, x);

pass(4) = ~isempty(F{1}) && ...
    norm(F1_vals - g_vals, inf) < 1e1*get(g, 'epslevel')*get(g, 'vscale');
pass(5) = ~isempty(F{2}) && ...
    norm(F2_vals - h_vals, inf) < max(get(h, 'epslevel')*get(h, 'vscale'));
pass(6) = ~isempty(F{3}) && ...
    norm(F3_vals - l_vals, inf) < get(l, 'epslevel')*get(l, 'vscale');

%% 
% [TODO]: Run a few tests for UNBNDFUN.
end
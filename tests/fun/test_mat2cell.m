% Test file for fun/mat2cell.m

function pass = test_mat2cell(pref)

% Get preferences.
if ( nargin < 2 )
    pref = chebpref();
end

% Set a domain for BNDFUN.
dom = [-2 7];

%% 
% Run a few tests for BNDFUN.

f = bndfun(@(x) [sin(x) cos(x) exp(x) x], dom, [], [], pref);
g = bndfun(@(x) sin(x), dom, [], [], pref);
h = bndfun(@(x) [cos(x) exp(x)], dom, [], [], pref);
l = bndfun(@(x) x, dom, [], [], pref);
    
% Test full arguments.
F = mat2cell(f, 1, [1 2 1]);
pass(n, 1) = ~isempty(F{1}) && normest(F{1} - g) < get(g, 'epslevel')*get(g, 'vscale');
pass(n, 2) = ~isempty(F{2}) && normest(F{2} - h) < max(get(h, 'epslevel')*get(h, 'vscale'));
pass(n, 3) = ~isempty(F{3}) && normest(F{3} - l) < get(l, 'epslevel')*get(l, 'vscale');
    
% Test two arguments.
F = mat2cell(f, [1 2 1]);
pass(n, 4) = ~isempty(F{1}) && normest(F{1} - g) < get(g, 'epslevel')*get(g, 'vscale');
pass(n, 5) = ~isempty(F{2}) && normest(F{2} - h) < max(get(h, 'epslevel')*get(h, 'vscale'));
pass(n, 6) = ~isempty(F{3}) && normest(F{3} - l) < get(l, 'epslevel')*get(l, 'vscale');
    
%% 
% [TODO]: Run a few tests for UNBNDFUN.
end
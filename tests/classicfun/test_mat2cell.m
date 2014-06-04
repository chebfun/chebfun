% Test file for @classicfun/mat2cell.m

function pass = test_mat2cell(pref)

% Get preferences.
if ( nargin < 2 )
    pref = chebfunpref();
end

% Set a domain for BNDFUN.
data.domain = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(data.domain) * rand(1000, 1) + data.domain(1);

%% 
% Run a few tests for BNDFUN.

f = bndfun(@(x) [sin(x) cos(x) exp(x) x], data, pref);
g = bndfun(@(x) sin(x), data, pref);
h = bndfun(@(x) [cos(x) exp(x)], data, pref);
l = bndfun(@(x) x, data, pref);

g_vals = feval(g, x);
h_vals = feval(h, x);
l_vals = feval(l, x);

% Test full arguments.
F = mat2cell(f, 1, [1 2 1]);
F1_vals = feval(F{1}, x);
F2_vals = feval(F{2}, x);
F3_vals = feval(F{3}, x);

err1 = normest(F{1} - g);
tol1 = 1e2*get(g, 'epslevel')*get(g, 'vscale');
pass(1) = ~isempty(F{1}) && (err1 < tol1);

err2 = normest(F{2} - h);
tol2 = 10*max(get(h, 'epslevel').*get(h, 'vscale'));
pass(2) = ~isempty(F{2}) && (err2 < tol2);

err3 = normest(F{3} - l);
tol3 = 10*get(l, 'epslevel')*get(l, 'vscale');
pass(3) = ~isempty(F{3}) && (err3 < tol3);

% Test two arguments.
F = mat2cell(f, [1 2 1]);
F1_vals = feval(F{1}, x);
F2_vals = feval(F{2}, x);
F3_vals = feval(F{3}, x);

err4 = normest(F{1} - g);
tol4 = 1e2*get(g, 'epslevel')*get(g, 'vscale');
pass(4) = ~isempty(F{1}) && (err4 < tol4);
    norm(F1_vals - g_vals, inf) < 1e1*get(g, 'epslevel')*get(g, 'vscale');

err5 = normest(F{2} - h);
tol5 = 10*max(get(h, 'epslevel').*get(h, 'vscale'));
pass(5) = ~isempty(F{2}) && (err5 < tol5);
    norm(F2_vals - h_vals, inf) < max(get(h, 'epslevel').*get(h, 'vscale'));

err6 = normest(F{3} - l);
tol6 = 10*get(l, 'epslevel')*get(l, 'vscale');
pass(6) = ~isempty(F{3}) && (err6 < tol6);

%% Test for UNBNDFUN:

% Functions on [-inf b]:

% Set the domain:
data.domain = [-Inf -3*pi];
domCheck = [-1e6 -3*pi];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

% Array-valued function:
op = @(x) [exp(x) x.*exp(x) (1-exp(x))./x];
opg = @(x) exp(x);
oph = @(x) [x.*exp(x) (1-exp(x))./x];

f = unbndfun(op, data);
F = mat2cell(f, 1, [1 2]);
F1Vals = feval(F{1}, x);
F2Vals = feval(F{2}, x);

F1Exact = opg(x);
F2Exact = oph(x);
err1 = F1Vals - F1Exact;
err2 = F2Vals - F2Exact;

pass(7) = norm([err1; err2(:)], inf) < ...
    1e1*max(get(f,'epslevel').*get(f,'vscale'));

end

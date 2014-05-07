function pass = test_extractColumns(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

f = chebfun(@(x) [sin(x), cos(x), exp(x)], pref);
g = chebfun(@(x) [sin(x), cos(x)], pref);
h = extractColumns(f, 1:2);
pass(1) = size(h, 2) == 2 && normest(g - h) < epslevel(f) && ...
    norm(g.pointValues - h.pointValues, inf) < epslevel(f);
h = f(:,1:2);
pass(2) = size(h, 2) == 2 && normest(g - h) < epslevel(f) && ...
    norm(g.pointValues - h.pointValues, inf) < epslevel(f);

g = chebfun(@(x) [sin(x), sin(x), exp(x), cos(x)], pref);
h = extractColumns(f, [1 1 3 2]);
pass(3) = size(h, 2) == 4 && normest(g - h) < epslevel(f) && ...
    norm(g.pointValues - h.pointValues, inf) < epslevel(f);

%% Test on function defined on unbounded domain:

% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -3*pi];
domCheck = [-1e6 -3*pi];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

% Array-valued function:
op = @(x) [exp(x) x.*exp(x) (1-exp(x))./x];
opg = @(x) [x.*exp(x) (1-exp(x))./x (1-exp(x))./x exp(x)];

f = chebfun(op, dom);
g = extractColumns(f, [2 3 3 1]);
gVals = feval(g, x);
gExact = opg(x);
err = gVals - gExact;

pass(4) = norm(err, inf) < 5*epslevel(f).*vscale(f);

end

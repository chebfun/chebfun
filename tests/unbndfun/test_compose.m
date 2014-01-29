% Test file for unbndfun/compose.m

function pass = test_compose(pref)

if ( nargin < 1 )
    pref = chebpref();
end

%%%%%%%%%%%%%%%%%% Compose an UNBNDFUN with an operator: %%%%%%%%%%%%%%%%%%

% Set the domain:
dom = [0 Inf];
domCheck = [0 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

opf = @(x) exp(-x);
opg = @(x) sin(exp(-x));
f = unbndfun(opf, dom);
g = compose(f, @sin);
gVals = feval(g, x);
gExact = opg(x);
err = gVals - gExact;
pass(1) = norm(err, inf) < get(g,'epslevel')*get(g,'vscale');







end
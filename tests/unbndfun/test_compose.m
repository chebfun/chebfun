% Test file for unbndfun/compose.m

function pass = test_compose(pref)

if ( nargin < 1 )
    pref = chebpref();
end

% Set the domain:
dom = [0 Inf];
domCheck = [0 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

%%%%%%%%%%%%%%%% Compose an UNBNDFUN with an operator (OP(F)) %%%%%%%%%%%%%%%%%%

opf = @(x) exp(-x);
opg = @(x) sin(exp(-x));
f = unbndfun(opf, dom);
g = compose(f, @sin);
gVals = feval(g, x);
gExact = opg(x);
err = gVals - gExact;
pass(1) = norm(err, inf) < get(g,'epslevel')*get(g,'vscale');

%%%%%%%%%%%%%%%%%%% Compose two UNBNDFUNs (F + G) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opf = @(x) exp(-x);
opg = @(x) x.*exp(-x);
oph = @(x) (x+1).*exp(-x);
f = unbndfun(opf, dom);
g = unbndfun(opg, dom);
h = compose(f, @plus, g);
hVals = feval(h, x);
hExact = oph(x);
err = hVals - hExact;
pass(2) = norm(err, inf) < get(h,'epslevel')*get(h,'vscale');

%%%%%%%%%%%%%%%%%%% Compose an UNBNDFUN with a BNDFUN (G(F)) %%%%%%%%%%%%%%%%%%%

opf = @(x) exp(-x);
opg = @(x) cos(x);
oph = @(x) cos(exp(-x));
f = unbndfun(opf, dom);
g = bndfun(opg, [-1 1]);
h = compose(f, g);
hVals = feval(h, x);
hExact = oph(x);
err = hVals - hExact;
pass(3) = norm(err, inf) < get(h,'epslevel')*get(h,'vscale');

end
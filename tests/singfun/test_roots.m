% Test file for singfun/roots.m

function pass = test_roots(pref)

% Get preferences.
if ( nargin < 1 )
    pref = singfun.pref;
end

% Set a tolerance.
tol = 10*pref.singfun.eps;

% The order of the exponents:
a = 0.56;
b = -0.56;
c = 1.28;
d = -1.28;

% Pre-allocate pass matrix
pass = zeros(1, 6);

% fractional root at the left endpoint and the smooth part has no roots in 
% [-1 1].
f = singfun(@(x) (1+x).^a.*exp(x), [a 0], {'root', 'none'}, pref);
r = roots(f);
r_exact = -1;
err = r - r_exact;
pass(1) = (norm(err, inf) < tol*normest(f));

% fractional pole at the left endpoint and note that the left endpoint is 
% not a root, even though the smooth part has a root there.
f = singfun(@(x) (1+x).^d.*sin(50*pi*x), [d+1 0], {'sing', 'none'}, pref);
r = roots(f);
r_exact = (-1+1/50:(1/50):1).';
err = r - r_exact;
pass(2) = (norm(err, inf) < tol*normest(f));

% fractional root at the right endpoint and the smooth part has no roots in [-1 1].
f = singfun(@(x) (1-x).^c.*cos(x), [0 c], {'none', 'sing'}, pref);
r = roots(f);
r_exact = 1;
err = r - r_exact;
pass(3) = (norm(err, inf) < tol*normest(f));

% no fractional pole but a root at the right endpoint.
f = singfun(@(x) (1-x).^b.*(exp(x)-exp(1)), [0 1+b], {'none', 'root'}, pref);
r = roots(f);
r_exact = 1;
err = r - r_exact;
pass(4) = (norm(err, inf) < tol*normest(f));

% a combination of fractional pole and fractional root.
f = singfun(@(x) (1+x).^b.*sin(x).*(1-x).^c, [b c], {'sing', 'root'}, pref);
r = roots(f);
r_exact = [0; 1];
err = r - r_exact;
pass(5) = (norm(err, inf) < tol*normest(f));

% Check the case with roots close to endpoints.
p = 1-1e-14;
f = singfun(@(x) (1+x).^b.*sin(x-p).*(1-x).^b, [b b], {'sing', 'sing'}, pref);
r = roots(f);
r_exact = p;
err = r - r_exact;
pass(6) = (norm(err, inf) < tol*normest(f));

end
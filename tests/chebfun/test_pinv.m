% Test file for @chebfun/pinv.m.

function pass = test_pinv(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Check a few simple examples.
f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1], pref);
g = pinv(f);
pass(1) = normest(f*(g*f) - f) < 1e2*vscale(f)*eps;

pass(2) = normest((g*f)*g - g) < 1e2*vscale(f)*eps;

f = chebfun(@(x) [sin(x) cos(x) exp(1i*x)], [-1 -0.5 0 0.5 1], pref);
g = pinv(f);
pass(3) = normest(f*(g*f) - f) < 1e2*vscale(f)*eps;

pass(4) = normest((g*f)*g - g) < 10*vscale(f)*eps;

% Check error conditions.
try
    g = pinv(f.');
    pass(5) = false;
catch ME
    pass(5) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:pinv:row');
end

end

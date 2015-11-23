% Test file for @chebfun/orth.m.

function pass = test_orth(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Check some simple examples.
f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 0 1], pref);
Q = orth(f);
pass(1) = norm(Q'*Q - eye(3), 'fro') < 10*vscale(Q)*eps;
pass(2) = normest(Q*(Q\f) - f) < 10*vscale(f)*eps;

f = chebfun(@(x) [sin(x) cos(x) exp(1i*x)], [-1 0 1], pref);
Q = orth(f);
pass(3) = norm(Q'*Q - eye(2), 'fro') < 10*vscale(Q)*eps;
pass(4) = normest(Q*(Q\f) - f) < 1e2*vscale(f)*eps;

% Check error conditions.
try
    Q = orth(f.');
    pass(5) = false;
catch ME
    pass(5) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:orth:row');
end

end

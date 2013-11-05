% Test file for @chebfun/orth.m.

function pass = test_orth(pref)

if ( nargin < 1 )
    pref = chebpref();
end

% Check some simple examples.
f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 0 1], pref);
Q = orth(f);
pass(1) = norm(Q'*Q - eye(3), 'fro') < 10*vscale(Q)*epslevel(Q);
pass(2) = normest(Q*(Q\f) - f) < 10*vscale(f)*epslevel(f);

f = chebfun(@(x) [sin(x) cos(x) exp(1i*x)], [-1 0 1], pref);
Q = orth(f);
pass(3) = norm(Q'*Q - eye(2), 'fro') < 10*vscale(Q)*epslevel(Q);
pass(4) = normest(Q*(Q\f) - f) < 10*vscale(f)*epslevel(f);

% Check error conditions.
try
    Q = orth(f.');
    pass(5) = false;
catch ME
    pass(5) = strcmp(ME.identifier, 'CHEBFUN:orth:row');
end

end

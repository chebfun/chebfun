% Test file for @chebfun/null.m.

function pass = test_null(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Check some simple examples.
f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 0 1], pref);
Z = null(f);
pass(1) = isempty(Z);

f = chebfun(@(x) [sin(x) cos(x) exp(1i*x)], [-1 0 1], pref);
Z = null(f);
pass(2) = isequal(size(Z), [3 1]);
pass(3) = abs(Z'*Z - 1) < 10*vscale(f)*eps;
pass(4) = normest(f*Z) < 10*vscale(f)*eps;

% Check error conditions.
try
    Q = null(f.');
    pass(5) = false;
catch ME
    pass(5) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:null:row');
end

end

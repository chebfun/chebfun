% Test file for @chebfun/svd.m.

function pass = test_svd()

if ( nargin < 1 )
    pref = chebfunpref();
end

% Check a simple example.
f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1]);
[U, S, V] = svd(f);
g = U*S*V';
pass(1) = ~U.isTransposed && (normest(f - g) < ...
    10*max(vscale(f)*epslevel(f), vscale(g)*epslevel(g)));
pass(2) = norm(U'*U - eye(3), 'fro') < 10*vscale(U)*epslevel(U);
pass(3) = norm(V'*V - eye(3), 'fro') < 10*vscale(f)*epslevel(f);
pass(4) = isequal(svd(f), diag(S));

ft = f.';
[U, S, V] = svd(ft);
g = U*S*V';
pass(4) = ~V.isTransposed && (normest(ft - g) < ...
    10*max(vscale(ft)*epslevel(ft), vscale(g)*epslevel(g)));
pass(5) = norm(U'*U - eye(3), 'fro') < 10*vscale(f)*epslevel(f);
pass(6) = norm(V'*V - eye(3), 'fro') < 10*vscale(V)*epslevel(V);

% Check error conditions.
try
    g = svd(f, 1);
    pass(7) = false;
catch ME
    pass(7) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:svd:twoArgs');
end

end

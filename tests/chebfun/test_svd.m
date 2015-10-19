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
    1e2*max(vscale(f)*eps, vscale(g)*eps));

pass(2) = norm(U'*U - eye(3), 'fro') < 10*vscale(U)*eps;
pass(3) = norm(V'*V - eye(3), 'fro') < 10*vscale(f)*eps;
pass(4) = isequal(svd(f), diag(S));

ft = f.';
[U, S, V] = svd(ft);
g = U*S*V';
pass(5) = ~V.isTransposed && (normest(ft - g) < ...
    1e2*max(vscale(ft)*eps, vscale(g)*eps));

pass(6) = norm(U'*U - eye(3), 'fro') < 10*vscale(f)*eps;
pass(7) = norm(V'*V - eye(3), 'fro') < 10*vscale(V)*eps;

% Check error conditions.
try
    g = svd(f, 1);
    pass(8) = false;
catch ME
    pass(7) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:svd:twoArgs');
end

end

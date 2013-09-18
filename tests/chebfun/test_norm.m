% Test file for @chebfun/norm.m.

function pass = test_norm(pref)

% Obtain preferences.
if ( nargin == 0 )
    pref = chebfun.pref();
end

% Check empty case.
pass(1) = norm(chebfun()) == 0;

% Check norms for scalar-valued chebfuns.
f = chebfun({@(x) exp(4*pi*1i*x), @exp, @exp}, [-1 0 0.5 1], pref);
normf_2 = sqrt(1 + (exp(2) - 1)/2);
pass(2) = abs(norm(f) - normf_2) < 10*vscale(f)*epslevel(f);
pass(3) = abs(norm(f, 2) - normf_2) < 10*vscale(f)*epslevel(f);
pass(4) = abs(norm(f, 'fro') - normf_2) < 10*vscale(f)*epslevel(f);

% [TODO]:  Test 1-norm.  (Needs abs().)

[maxVal, maxLoc] = norm(f, inf);
pass(5) = abs(maxVal - exp(1)) < 10*vscale(f)*epslevel(f) && ...
    abs(feval(f, maxLoc) - exp(1)) < 10*vscale(f)*epslevel(f);

[minVal, minLoc] = norm(f, -inf);
pass(6) = abs(minVal - 1) < 10*vscale(f)*epslevel(f) && ...
    abs(feval(f, minLoc) - 1) < 10*vscale(f)*epslevel(f);

g = chebfun(@(x) 1./(1 + (x - 0.1).^2), [-1 -0.5 0 0.5 1]);
[maxVal, maxLoc] = norm(g, inf);
pass(7) = abs(maxVal - 1) < 10*vscale(g)*epslevel(g) && ...
    abs(feval(g, maxLoc) - 1) < 10*vscale(g)*epslevel(g);

% [TODO]:  Test p-norm.  (Needs power().)

% Check norms for array-valued chebfuns.
f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1], pref);

% [TODO]:  Test 1-norm.  (Needs abs().)

% [TODO]:  Test 2-norm.  (Needs svd().)

pass(8) = abs(norm(f) - 2.372100421113536830) < 10*vscale(f)*epslevel(f);
pass(9) = abs(norm(f, 'fro') - 2.372100421113536830) < 10*vscale(f)*epslevel(f);

% [TODO]:  Test inf norm.  (Needs abs().)

% [TODO]:  Test -inf "norm".  (Needs abs().)

% [TODO]:  Test p-norm.  (Needs power().)

% Check error conditions.
try
    [x, y] = norm(g, 1);
    pass(10) = false;
catch ME
    pass(10) = strcmp(ME.identifier, 'CHEBFUN:norm:argout');
end

try
    [x, y] = norm(g, 2);
    pass(11) = false;
catch ME
    pass(11) = strcmp(ME.identifier, 'CHEBFUN:norm:argout');
end

try
    [x, y] = norm(g, 0.4);
    pass(12) = false;
catch ME
    pass(12) = strcmp(ME.identifier, 'CHEBFUN:norm:argout');
end

try
    [x, y] = norm(g, 'bad');
    pass(13) = false;
catch ME
    pass(13) = strcmp(ME.identifier, 'CHEBFUN:norm:unknownNorm');
end

try
    [x, y] = norm(f, 2);
    pass(14) = false;
catch ME
    pass(14) = strcmp(ME.identifier, 'CHEBFUN:norm:argout');
end

try
    [x, y] = norm(f, 'fro');
    pass(15) = false;
catch ME
    pass(15) = strcmp(ME.identifier, 'CHEBFUN:norm:argout');
end

try
    [x, y] = norm(f, 'bad');
    pass(16) = false;
catch ME
    pass(16) = strcmp(ME.identifier, 'CHEBFUN:norm:unknownNorm');
end

end

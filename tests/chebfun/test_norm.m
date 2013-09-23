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
pass(5) = abs(norm(f, 1) - exp(1)) < 10*vscale(f)*epslevel(f);

[maxVal, maxLoc] = norm(f, inf);
pass(6) = abs(maxVal - exp(1)) < 10*vscale(f)*epslevel(f) && ...
    abs(feval(f, maxLoc) - exp(1)) < 10*vscale(f)*epslevel(f);

[minVal, minLoc] = norm(f, -inf);
pass(7) = abs(minVal - 1) < 10*vscale(f)*epslevel(f) && ...
    abs(feval(f, minLoc) - 1) < 10*vscale(f)*epslevel(f);

g = chebfun(@(x) 1./(1 + (x - 0.1).^2), [-1 -0.5 0 0.5 1]);
[maxVal, maxLoc] = norm(g, inf);
pass(8) = abs(maxVal - 1) < 10*vscale(g)*epslevel(g) && ...
    abs(feval(g, maxLoc) - 1) < 10*vscale(g)*epslevel(g);

% [TODO]:  Test p-norm.  (Needs power().)

% Check norms for array-valued chebfuns.
f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1], pref);

[normVal, col] = norm(f, 1);
pass(9) = (col == 3) && abs(normVal - (exp(1) - exp(-1))) ...
    < 10*vscale(f)*epslevel(f);

% [TODO]:  Test 2-norm.  (Needs svd().)

pass(10) = abs(norm(f) - 2.372100421113536830) < 10*vscale(f)*epslevel(f);
pass(11) = abs(norm(f, 'fro') - 2.372100421113536830) < ...
    10*vscale(f)*epslevel(f);

[normVal, loc] = norm(f, inf);
pass(12) = (loc == 1) && abs(normVal - (exp(1) + sin(1) + cos(1))) ...
    < 10*vscale(f)*epslevel(f);

[normVal, loc] = norm(f, -inf);
pass(13) = (loc == -1) && abs(normVal - (exp(-1) + sin(1) + cos(1))) ...
    < 10*vscale(f)*epslevel(f);

% [TODO]:  Test p-norm.  (Needs power().)

% Check error conditions.
try
    [x, y] = norm(g, 1);
    pass(14) = false;
catch ME
    pass(14) = strcmp(ME.identifier, 'CHEBFUN:norm:argout');
end

try
    [x, y] = norm(g, 2);
    pass(15) = false;
catch ME
    pass(15) = strcmp(ME.identifier, 'CHEBFUN:norm:argout');
end

try
    [x, y] = norm(g, 0.4);
    pass(16) = false;
catch ME
    pass(16) = strcmp(ME.identifier, 'CHEBFUN:norm:argout');
end

try
    [x, y] = norm(g, 'bad');
    pass(17) = false;
catch ME
    pass(17) = strcmp(ME.identifier, 'CHEBFUN:norm:unknownNorm');
end

try
    [x, y] = norm(f, 2);
    pass(18) = false;
catch ME
    pass(18) = strcmp(ME.identifier, 'CHEBFUN:norm:argout');
end

try
    [x, y] = norm(f, 'fro');
    pass(19) = false;
catch ME
    pass(19) = strcmp(ME.identifier, 'CHEBFUN:norm:argout');
end

try
    [x, y] = norm(f, 'bad');
    pass(20) = false;
catch ME
    pass(20) = strcmp(ME.identifier, 'CHEBFUN:norm:unknownNorm');
end

end

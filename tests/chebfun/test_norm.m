% Test file for @chebfun/norm.m.

function pass = test_norm(pref)

% Obtain preferences.
if ( nargin == 0 )
    pref = chebpref();
end

% Check empty case.
pass(1) = norm(chebfun()) == 0;

% Check norms for scalar-valued chebfuns.
f = chebfun({@(x) exp(4*pi*1i*x), @exp, @exp}, [-1 0 0.5 1], pref);
normf_2 = sqrt(1 + (exp(2) - 1)/2);
pass(2) = abs(norm(f) - normf_2) < 10*vscale(f).*epslevel(f);
pass(3) = abs(norm(f, 2) - normf_2) < 10*vscale(f).*epslevel(f);
pass(4) = abs(norm(f, 'fro') - normf_2) < 10*vscale(f).*epslevel(f);
pass(5) = abs(norm(f, 1) - exp(1)) < 10*vscale(f).*epslevel(f);

[maxVal, maxLoc] = norm(f, inf);
pass(6) = abs(maxVal - exp(1)) < 10*vscale(f).*epslevel(f) && ...
    abs(feval(f, maxLoc) - exp(1)) < 10*vscale(f).*epslevel(f);

[minVal, minLoc] = norm(f, -inf);
pass(7) = abs(minVal - 1) < 10*vscale(f).*epslevel(f) && ...
    abs(feval(f, minLoc) - 1) < 10*vscale(f).*epslevel(f);

g = chebfun(@(x) 1./(1 + (x - 0.1).^2), [-1 -0.5 0 0.5 1]);
[maxVal, maxLoc] = norm(g, inf);
pass(8) = abs(maxVal - 1) < 10*vscale(g).*epslevel(g) && ...
    abs(feval(g, maxLoc) - 1) < 10*vscale(g).*epslevel(g);

% Check norms for array-valued chebfuns.
f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1], pref);

[normVal, col] = norm(f, 1);
pass(9) = (col == 3) && abs(normVal - (exp(1) - exp(-1))) ...
    < 10*vscale(f).*epslevel(f);

U = chebfun(@(x) [(1 + 0*x) exp(2*pi*1i*x) exp(2*pi*1i*2*x)], [0 1], pref);
S = diag([pi ; exp(1) ; 1]);
V = [1/sqrt(2) -1/sqrt(2) 0 ; 1/sqrt(2) 1/sqrt(2) 0 ; 0 0 1];
h = U*S*V';
pass(10) = abs(norm(h, 2) - pi) < 10*vscale(h).*epslevel(h);

pass(11) = abs(norm(f) - 2.372100421113536830) < 10*vscale(f).*epslevel(f);
pass(12) = abs(norm(f, 'fro') - 2.372100421113536830) < ...
    10*vscale(f).*epslevel(f);

[normVal, loc] = norm(f, inf);
pass(13) = (loc == 1) && abs(normVal - (exp(1) + sin(1) + cos(1))) ...
    < 10*vscale(f).*epslevel(f);

[normVal, loc] = norm(f, -inf);
pass(14) = (loc == -1) && abs(normVal - (exp(-1) + sin(1) + cos(1))) ...
    < 10*vscale(f).*epslevel(f);

% Check error conditions.
try
    [x, y] = norm(g, 1);
    pass(15) = false;
catch ME
    pass(15) = strcmp(ME.identifier, 'CHEBFUN:norm:argout');
end

try
    [x, y] = norm(g, 2);
    pass(16) = false;
catch ME
    pass(16) = strcmp(ME.identifier, 'CHEBFUN:norm:argout');
end

try
    [x, y] = norm(g, 0.4);
    pass(17) = false;
catch ME
    pass(17) = strcmp(ME.identifier, 'CHEBFUN:norm:argout');
end

try
    [x, y] = norm(g, 'bad');
    pass(18) = false;
catch ME
    pass(18) = strcmp(ME.identifier, 'CHEBFUN:norm:unknownNorm');
end

try
    [x, y] = norm(f, 2);
    pass(19) = false;
catch ME
    pass(19) = strcmp(ME.identifier, 'CHEBFUN:norm:argout');
end

try
    [x, y] = norm(f, 'fro');
    pass(20) = false;
catch ME
    pass(20) = strcmp(ME.identifier, 'CHEBFUN:norm:argout');
end

try
    [x, y] = norm(f, 'bad');
    pass(21) = false;
catch ME
    pass(21) = strcmp(ME.identifier, 'CHEBFUN:norm:unknownNorm');
end

%% Tests on SINGFUN-related stuff:

% p-norm (p = 3) - Tests on power:
f = chebfun(@(x) sin(50*x), 'splitting', 'on');
p = norm(f, 3);
p_exact = 0.9484869030456855; % This is obtained using Mathematica.
err = p-p_exact;
pass(22) = norm(err, inf) < vscale(f).*epslevel(f);

% Test on finite SINGFUN:

% 2-norm:
f = chebfun(@(x) sin(100*x).*((x+1).^0.6), 'exps', [0.6 0], 'splitting', 'on');
p2 = norm(f, 2);
p2_exact = 1.02434346249849423; % This is obtained using Mathematica.
err = p2-p2_exact;
pass(23) = norm(err, inf) < vscale(f).*epslevel(f);

% 1-norm:
p1 = norm(f, 1);
p1_exact = 1.20927413792339491; % This is obtained using Mathematica.
err = p1-p1_exact;
pass(24) = norm(err, inf) < 1e5*vscale(f).*epslevel(f);

% Inf-norm:
f = chebfun(@(x) sin(x).*((1-x).^0.6), 'exps', [0 0.6], 'splitting', 'on');
[normF, normLoc] = norm(f, Inf);
p_exact = [1.275431511911148 -1];
err = [normF, normLoc] - p_exact;
pass(25) = norm(err, inf) < vscale(f).*epslevel(f);

% -Inf-norm:
f = chebfun(@(x) (sin(x)-0.4).*((x+1).^0.8), 'exps', [0.8 0], 'splitting', 'on');
[normF, normLoc] = norm(f, -Inf);
p_exact = [0 -1];
err = [normF, normLoc] - p_exact;
pass(26) = norm(err, inf) < vscale(f).*epslevel(f);

end
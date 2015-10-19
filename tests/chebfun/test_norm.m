% Test file for @chebfun/norm.m.

function pass = test_norm(pref)

% Obtain preferences.
if ( nargin == 0 )
    pref = chebfunpref();
end

% Check empty case.
pass(1) = norm(chebfun()) == 0;

% Check norms for scalar-valued chebfuns.
f = chebfun({@(x) exp(4*pi*1i*x), @exp, @exp}, [-1 0 0.5 1], pref);
normf_2 = sqrt(1 + (exp(2) - 1)/2);
pass(2) = abs(norm(f) - normf_2) < 10*vscale(f)*eps;
pass(3) = abs(norm(f, 2) - normf_2) < 10*vscale(f)*eps;
pass(4) = abs(norm(f, 'fro') - normf_2) < 10*vscale(f)*eps;
pass(5) = abs(norm(f, 1) - exp(1)) < 10*vscale(f)*eps;

[maxVal, maxLoc] = norm(f, inf);
pass(6) = abs(maxVal - exp(1)) < 10*vscale(f)*eps && ...
    abs(feval(f, maxLoc) - exp(1)) < 10*vscale(f)*eps;

[minVal, minLoc] = norm(f, -inf);
pass(7) = abs(minVal - 1) < 10*vscale(f)*eps && ...
    abs(abs(feval(f, minLoc)) - 1) < 10*vscale(f)*eps;

g = chebfun(@(x) 1./(1 + (x - 0.1).^2), [-1 -0.5 0 0.5 1]);
[maxVal, maxLoc] = norm(g, inf);
pass(8) = abs(maxVal - 1) < 10*vscale(g)*eps && ...
    abs(feval(g, maxLoc) - 1) < 10*vscale(g)*eps;

% Check norms for array-valued chebfuns.
f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1], pref);

[normVal, col] = norm(f, 1);
pass(9) = (col == 3) && abs(normVal - (exp(1) - exp(-1))) ...
    < 10*vscale(f)*eps;

U = chebfun(@(x) [(1 + 0*x) exp(2*pi*1i*x) exp(2*pi*1i*2*x)], [0 1], pref);
S = diag([pi ; exp(1) ; 1]);
V = [1/sqrt(2) -1/sqrt(2) 0 ; 1/sqrt(2) 1/sqrt(2) 0 ; 0 0 1];
h = U*S*V';
pass(10) = abs(norm(h, 2) - pi) < 10*vscale(h)*eps;

pass(11) = abs(norm(f) - 2.372100421113536830) < 10*vscale(f)*eps;
pass(12) = abs(norm(f, 'fro') - 2.372100421113536830) < ...
    10*vscale(f)*eps;

[normVal, loc] = norm(f, inf);
pass(13) = (loc == 1) && abs(normVal - (exp(1) + sin(1) + cos(1))) ...
    < 10*vscale(f)*eps;

[normVal, loc] = norm(f, -inf);
pass(14) = (loc == -1) && abs(normVal - (exp(-1) + sin(1) + cos(1))) ...
    < 10*vscale(f)*eps;

% Check error conditions.
try
    [x, y] = norm(g, 1);
    pass(15) = false;
catch ME
    pass(15) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:norm:argout');
end

try
    [x, y] = norm(g, 2);
    pass(16) = false;
catch ME
    pass(16) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:norm:argout');
end

try
    [x, y] = norm(g, 0.4);
    pass(17) = false;
catch ME
    pass(17) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:norm:argout');
end

try
    [x, y] = norm(g, 'bad');
    pass(18) = false;
catch ME
    pass(18) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:norm:unknownNorm');
end

try
    [x, y] = norm(f, 2);
    pass(19) = false;
catch ME
    pass(19) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:norm:argout');
end

try
    [x, y] = norm(f, 'fro');
    pass(20) = false;
catch ME
    pass(20) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:norm:argout');
end

try
    [x, y] = norm(f, 'bad');
    pass(21) = false;
catch ME
    pass(21) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:norm:unknownNorm');
end

%% Tests for singular functions:

% p-norm (p = 3) - Tests on power:
f = chebfun(@(x) sin(50*x), 'splitting', 'on');
p = norm(f, 3);
p_exact = 0.9484869030456855; % This is obtained using Mathematica.
err = p-p_exact;
pass(22) = norm(err, inf) < 1e1*vscale(f)*eps;


% Test on finite SINGFUN:

% 2-norm:
f = chebfun(@(x) sin(100*x).*((x+1).^0.6), 'exps', [0.6 0], 'splitting', 'on');
p2 = norm(f, 2);
p2_exact = 1.02434346249849423; % This is obtained using Mathematica.
err = p2-p2_exact;
pass(23) = norm(err, inf) < 1e1*vscale(f)*eps;


% 1-norm:
p1 = norm(f, 1);
p1_exact = 1.20927413792339491; % This is obtained using Mathematica.
err = p1-p1_exact;
pass(24) = norm(err, inf) < 1e7*vscale(f)*eps;


% Inf-norm:
f = chebfun(@(x) sin(x).*((1-x).^0.6), 'exps', [0 0.6], 'splitting', 'on');
[normF, normLoc] = norm(f, Inf);
p_exact = [1.275431511911148 -1]; % This is obtained using Mathematica.
err = [normF, normLoc] - p_exact;
pass(25) = norm(err, inf) < 1e1*vscale(f)*eps;


% -Inf-norm:
f = chebfun(@(x) (sin(x+1.1)).*((x+1).^0.8), 'exps', [0.8 0], 'splitting', 'on');
[normF, normLoc] = norm(f, -Inf);
p_exact = [0 -1]; % This is obtained using Mathematica.
err = [normF, normLoc] - p_exact;
pass(26) = norm(err, inf) < 1e3*vscale(f)*eps;
%% Test for functions defined on unbounded domain:

% Functions on [a inf]:

% Set the domain:
dom = [1 Inf];

% 2-norm and fro-norm:
op = @(x) (1-exp(-x))./x;
f = chebfun(op, dom);
p = norm(f);
pExact = 0.860548225417173;  % This is obtained using Matlab symbolic toolbox.
err = p - pExact;
pass(27) = abs(err) < 1e7*eps*vscale(f);

% 1-norm:
op = @(x) (1-exp(-x))./(x.^2);
f = chebfun(op, dom);
p = norm(f, 1);
pExact = 0.851504493224078;  % This is obtained using Matlab symbolic toolbox.
err = p - pExact;
pass(28) = abs(err) < 1e6*eps*vscale(f);

% P-norm (here P = 3):
op = @(x) (1-exp(-x))./x;
f = chebfun(op, dom);
p = norm(f, 3);
pExact = 0.631964633132246;  % This is obtained using Matlab symbolic toolbox.
err = p - pExact;
pass(29) = abs(err) < 1e6*eps*vscale(f);

% Inf-norm:
[normF, normLoc] = norm(f, Inf);
p_exact = [op(1) 1];
err = [normF, normLoc] - p_exact;
pass(30) = norm(err, inf) < vscale(f)*eps;

% -Inf-norm:
op = @(x) 1./x;
f = chebfun(op, dom);
[normF, normLoc] = norm(f, -Inf);
pass(31) = ( abs(normF) < vscale(f)*eps ) && ...
    ( normLoc == Inf );

%% Array-valued function defined on unbounded domain:

% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -1];

% 1-norm:
op = @(x) [exp(x) x.*exp(x) (1-exp(x))./(x.^2)];

f = chebfun(op, dom);
p = norm(f{-1000,-1}, 1);
pExact = 0.851504493224078;  % This is obtained using Matlab symbolic toolbox.
err = p - pExact;
% The tolerance below is loosen to allow certain spurious roots:
pass(32) = abs(err) < 1e-2;

%% #578
f = chebfun(@(x) cos(x)./(1e5+(x-30).^6),[0 inf]);
warning('off', 'CHEBFUN:UNBNDFUN:sum:slowdecay')
I = norm(f);
warning('on', 'CHEBFUN:UNBNDFUN:sum:slowdecay')
% The following result is obtained using Mathematica:
Iexact = 2.4419616835794597e-5;
err = abs(I - Iexact);
pass(33) = ( err < 1e9*vscale(f)*eps );

% #920: sum of array-valued chebfun defined on unbounded domain:
% (same function which decays fast enough to be integrable):
f = chebfun(@(x) exp(-[x x].^2), [0 Inf]);
I = norm(f);
% The following exact value is obtained by Mathematica.
I_exact = 1.119515134920248;
pass(34) = norm(I-I_exact, inf) < 1e1*eps*get(f,'vscale');

end

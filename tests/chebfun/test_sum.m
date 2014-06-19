% Test file for @chebfun/sum.m.

function pass = test_sum(pref)

% Obtain preferences.
if ( nargin == 0 )
    pref = chebfunpref();
end

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

%% SCALAR-VALUED

% Check the empty case.
pass(1) = sum(chebfun()) == 0;

% Check operation in the general case.
f = chebfun({@(x) exp(4*pi*1i*x), @exp, @exp}, [-1 0 0.5 1], pref);
pass(2) = abs(sum(f) - (exp(1) - 1)) < 10*vscale(f)*epslevel(f);

% Check operation for row chebfuns.
ft = f.';
pass(3) = abs(sum(ft) - (exp(1) - 1)) < 10*vscale(ft)*epslevel(ft);

% Check sum over a subdomain.
pass(4) = abs(sum(f, [-1 1]) - (exp(1) - 1)) < 10*vscale(f)*epslevel(f);
pass(5) = abs(sum(f, [-1 0])) < 10*vscale(f)*epslevel(f);
pass(6) = abs(sum(f, [0 1]) - (exp(1) - 1)) < 10*vscale(f)*epslevel(f);

% Check sum between chebfun limits.
f = chebfun(@exp, [-1 -0.5 0 0.5 1], pref);
a = chebfun(@(x) x.^2 - 1, [-1 1]);
b = chebfun(@(x) -x.^2 + 1, [-1 1]);

F1 = sum(f, a, 1);
F1_exact = @(x) exp(1) - exp(x.^2 - 1);
pass(7) = norm(feval(F1, xr) - F1_exact(xr), inf) < 10*vscale(F1)*epslevel(F1);

F2 = sum(f, -1, b);
F2_exact = @(x) exp(-x.^2 + 1) - exp(-1);
pass(8) = norm(feval(F2, xr) - F2_exact(xr), inf) < 10*vscale(F2)*epslevel(F2);

F3 = sum(f, a, b);
F3_exact = @(x) exp(-x.^2 + 1) - exp(x.^2 - 1);
pass(9) = norm(feval(F3, xr) - F3_exact(xr), inf) < 10*vscale(F3)*epslevel(F3);

%% ARRAY-VALUED
f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1]);
pass(10) = norm(sum(f) - [0 2*sin(1) (exp(1) - exp(-1))], inf) < ...
    10*vscale(f)*epslevel(f);

ft = f.';
pass(11) = norm(sum(ft) - [0 2*sin(1) (exp(1) - exp(-1))], inf) < ...
    10*vscale(ft)*epslevel(ft);

pass(12) = norm(sum(f, [-1 1]) - [0 2*sin(1) (exp(1) - exp(-1))], inf) < ...
    10*vscale(f)*epslevel(f);

pass(13) = norm(sum(f, [-1 0]) - [(cos(-1) - 1) sin(1) (1 - exp(-1))], inf) ...
    < 10*vscale(f)*epslevel(f);

F1 = sum(f, a, 1);
F1_col1_exact = @(x) cos(x.^2 - 1) - cos(1);
F1_col2_exact = @(x) sin(1) - sin(x.^2 - 1);
F1_col3_exact = @(x) exp(1) - exp(x.^2 - 1);
F1_exact = @(x) [F1_col1_exact(x) F1_col2_exact(x) F1_col3_exact(x)];
err = feval(F1, xr) - F1_exact(xr);
pass(14) = norm(err(:), inf) < 10*vscale(F1)*epslevel(F1);

F2 = sum(f, -1, b);
F2_col1_exact = @(x) cos(-1) - cos(-x.^2 + 1);
F2_col2_exact = @(x) sin(-x.^2 + 1) - sin(-1);
F2_col3_exact = @(x) exp(-x.^2 + 1) - exp(-1);
F2_exact = @(x) [F2_col1_exact(x) F2_col2_exact(x) F2_col3_exact(x)];
err = feval(F2, xr) - F2_exact(xr);
pass(15) = norm(err(:), inf) < 10*vscale(F2)*epslevel(F2);

F3 = sum(f, a, b);
F3_col1_exact = @(x) -cos(-x.^2 + 1) + cos(x.^2 - 1);
F3_col2_exact = @(x) sin(-x.^2 + 1) - sin(x.^2 - 1);
F3_col3_exact = @(x) exp(-x.^2 + 1) - exp(x.^2 - 1);
F3_exact = @(x) [F3_col1_exact(x) F3_col2_exact(x) F3_col3_exact(x)];
err = feval(F3, xr) - F3_exact(xr);
pass(16) = norm(err(:), inf) < 10*vscale(F3)*epslevel(F3);

%% Check dim argument.
g = sum(f, 2);
g_exact = @(x) sin(x) + cos(x) + exp(x);
pass(17) = norm(feval(g, xr) - g_exact(xr), inf) < 10*vscale(g)*epslevel(g);

g = sum(ft, 1);
g_exact = @(x) (sin(x) + cos(x) + exp(x)).';
pass(18) = norm(feval(g, xr) - g_exact(xr), inf) < 10*vscale(g)*epslevel(g);

%% Check error conditions.
try
    s = sum(f, -2, 2);
    pass(19) = false;
catch ME
    pass(19) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:sum:sumSubDom:ab');
end

try
    s = sum(f, -2, b);
    pass(20) = false;
catch ME
    pass(20) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:sum:sumSubDom:a');
end

try
    s = sum(f, a, 2);
    pass(21) = false;
catch ME
    pass(21) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:sum:sumSubDom:b');
end

%% QUASIMATRICES:

f = quasimatrix(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1]);
pass(22) = norm(sum(f) - [0 2*sin(1) (exp(1) - exp(-1))], inf) < ...
    10*vscale(f)*epslevel(f);

ft = f.';
pass(23) = norm(sum(ft) - [0 2*sin(1) (exp(1) - exp(-1))], inf) < ...
    10*vscale(ft)*epslevel(ft);

pass(24) = norm(sum(f, [-1 1]) - [0 2*sin(1) (exp(1) - exp(-1))], inf) < ...
    10*vscale(f)*epslevel(f);

pass(25) = norm(sum(f, [-1 0]) - [(cos(-1) - 1) sin(1) (1 - exp(-1))], inf) ...
    < 10*vscale(f)*epslevel(f);

%% Test on singular function: piecewise smooth chebfun - splitting on.

% Set a domain
dom = [-2 7];

pow = -0.5;
op = @(x) (x - dom(1)).^pow.*sin(100*x);
f = chebfun(op, dom, 'exps', [pow 0], 'splitting', 'on');
I = sum(f);
I_exact = 0.17330750941063138;
pass(26) = ( abs(I-I_exact) < 20*get(f, 'epslevel')*abs(I_exact) );

%% Test for functions defined on unbounded domain:

% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf 2 Inf];

op1 = @(x) x.^2.*exp(-x.^2);
op2 = @(x) (1-exp(-x.^2))./x.^2;
f = chebfun({op1 op2}, dom);
I = sum(f);
% The exact solution is obtained using Matlab symbolic toolbox:
IExact = 1.364971769155161;
err = abs(I - IExact);
pass(27) = err < 1e7*get(f,'epslevel')*get(f,'vscale');

% Functions on [0 inf]:

% Set the domain:
dom = [0 Inf];

op1 = @(x) x.*exp(-x);
f = chebfun(op1, dom);
I1 = sum(f);

IExact = 1;
err1 = abs(I1 - IExact);
pass(28) = err1 < 2e5*get(f,'epslevel')*get(f,'vscale');

x = chebfun('x', dom);
g = x.*exp(-x);
warning('off', 'CHEBFUN:UNBNDFUN:sum:slowDecay');
I2 = sum(g);
warning('on', 'CHEBFUN:UNBNDFUN:sum:slowdDecay');
err2 = abs(I2 - IExact);
tol = 2e10*get(f,'epslevel')*get(f,'vscale');
pass(29) = err2 < tol;

% Function on [-Inf Inf]:
f = chebfun('exp(-x.^2/16).*(1+.2*cos(10*x))',[-inf,inf]);

% Suppress expected warnings which may occur on certain machines:
warning('off','CHEBFUN:UNBNDFUN:sum:slowDecay'); 
I = sum(f);
% Re-enable warnings:
warning('on','CHEBFUN:UNBNDFUN:sum:slowDecay');  
IExact = 7.0898154036220641;
err = abs(I - IExact);
pass(30) = err < 1e9*get(f,'epslevel')*get(f,'vscale');

% #496
k = 3;
op = @(t) exp(-t)./(1+k*t);
f = chebfun(op, [0 inf]);
I = sum(f);
% The following exact result is obtained using Matlab symbolic toolbox:
I_exact = 0.385602012136694;
err = abs(I-I_exact);
pass(31) = err < 2e1*get(f,'epslevel')*get(f,'vscale');

% #920: sum of array-valued chebfun defined on unbounded domain:
% (same function which decays fast enough to be integrable):
f = chebfun(@(x) exp(-[x x].^2), [0 Inf]);
I = sum(f);
% The following exact value is obtained by Mathematica.
I_exact = 0.886226925452758*ones(1, 2);
pass(32) = norm(I-I_exact, inf) < 1e1*get(f,'epslevel')*get(f,'vscale');

% #920: sum of array-valued chebfun defined on unbounded domain: 
% (different functions but both decay fast so that integrable)
f = chebfun(@(x)[exp(-x.^2) 1./(x.^3)], [1 Inf]);
I = sum(f);
% The following exact value is obtained by Mathematica.
I_exact = [0.139402792640331 0.5];
pass(33) = norm(I-I_exact, inf) < 1e6*get(f,'epslevel')*get(f,'vscale');

% #920: sum of array-valued chebfun defined on unbounded domain: 
% (two different functions one of which is integrable and the other one is not):
f = chebfun(@(x)[exp(-x.^2) 0*x+1], [1 Inf]);
I = sum(f);
% The following exact value is obtained by Mathematica.
I_exact = 0.139402792640331;
pass(34) = abs(I(1)-I_exact) < get(f,'epslevel')*get(f,'vscale') && ...
    isinf(I(2));

end
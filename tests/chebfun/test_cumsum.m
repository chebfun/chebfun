% Test file for @chebfun/cumsum.m.

function pass = test_cumsum(pref)

% Obtain preferences.
if ( nargin == 0 )
    pref = chebfunpref();
end

%%
% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

% Check empty case.
pass(1) = isempty(cumsum(chebfun()));

%%
% Check an example with breakpoints but no impulses.
f1 = chebfun(@cos, [-1 -0.5 0.5 1], pref);
If1 = cumsum(f1);
If1_exact = @(x) sin(x) - sin(-1);
pass(2) = norm(feval(If1, xr) - If1_exact(xr), inf) < ...
    10*vscale(If1)*epslevel(If1);

%%
% Check behavior for row chebfuns.
f1t = f1.';
If1t = cumsum(f1t);
If1t_exact = @(x) (sin(x) - sin(-1)).';
pass(3) = norm(feval(If1t, xr) - If1t_exact(xr), inf) < ...
    10*vscale(If1t)*epslevel(If1t);

%%
% Check behavior for array-valued chebfuns.
f3 = chebfun(@(x) [cos(x) -sin(x) exp(x)], [-1 -0.5 0.5 1], pref);
If3 = cumsum(f3);
If3_exact = @(x) [sin(x) cos(x) exp(x)] - ...
    repmat([sin(-1) cos(-1) exp(-1)], size(x, 1), 1);
pass(4) = max(max(abs(feval(If3, xr) - If3_exact(xr)))) < ...
    10*vscale(If3)*epslevel(If3);

f3t = f3.';
If3t = cumsum(f3t);
If3t_exact = @(x) ([sin(x) cos(x) exp(x)] - ...
    repmat([sin(-1) cos(-1) exp(-1)], size(x, 1), 1)).';
pass(5) = max(max(abs(feval(If3t, xr) - If3t_exact(xr)))) < ...
    10*vscale(If3t)*epslevel(If3t);

f4 = chebfun(@(x) [cos(x) -sin(x)], [-1 -0.5 0.5 1], pref);
If4 = cumsum(f4);
If4_exact = @test_If4;
pass(6) = max(max(abs(feval(If4, xr) - If4_exact(xr)))) < ...
    10*vscale(If4)*epslevel(If4);

% Check second argument.
I2f1 = cumsum(f1, 2);
I2f1_exact = @(x) -cos(x) - sin(-1)*x - (-cos(-1) - sin(-1)*-1);
pass(7) = norm(feval(I2f1, xr) - I2f1_exact(xr), inf) < ...
    10*vscale(I2f1)*epslevel(I2f1);

% Check second argument.
I2f1 = cumsum(f1, 2);
I2f1_exact = @(x) -cos(x) - sin(-1)*x - (-cos(-1) - sin(-1)*-1);
pass(8) = norm(feval(I2f1, xr) - I2f1_exact(xr), inf) < ...
    10*vscale(I2f1)*epslevel(I2f1);

I2f3 = cumsum(f3, 2);
I2f3_exact = @(x) [-cos(x) sin(x) exp(x)] - ...
    [sin(-1)*x cos(-1)*x exp(-1)*x] - ...
    repmat([(-cos(-1) - sin(-1)*-1) (sin(-1) - cos(-1)*-1) ...
            (exp(-1) - exp(-1)*-1)], size(x, 1), 1);
pass(9) = max(max(abs(feval(I2f3, xr) - I2f3_exact(xr)))) < ...
    10*vscale(I2f3)*epslevel(I2f3);

%% Test on singular function:

dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

pow = -1.64;
op = @(x) (x-dom(1)).^pow;
f = chebfun(op, dom, 'exps', [pow 0]);
g = cumsum(f);
vals_g = feval(g, x); 
g_exact = @(x) (x-dom(1)).^(pow+1)./(pow+1);
vals_exact = feval(g_exact, x);
err = vals_g - vals_exact;
pass(10) = ( norm(err, inf) < 1e3*get(f,'epslevel')*norm(vals_exact, inf) );

%% piecewise smooth chebfun: smoothfun + singfun & splitting off:

% define the domain:
dom = [-2 -1 0 1];

op1 = @(x) sin(x);
op2 = @(x) 1./((1+x).^0.5);
op3 = @(x) x+1;
opi1 = @(x) cos(-2)-cos(x);
opi2 = @(x) 2.*((1+x).^0.5);
opi3 = @(x) x.^2/2 + x + opi2(0);

op = {op1, op2, op3};
opi = {opi1, opi2, opi3};
f = chebfun(op, dom, 'exps', [0 0 -0.5 0 0 0]);

% We temporarily disable this warning: 
warning('off', 'CHEBFUN:SINGFUN:plus:exponentDiff');
g = cumsum(f);
warning('on', 'CHEBFUN:SINGFUN:plus:exponentDiff');

% check values:
result = zeros(1,3);
for j = 1:3
    % define the domain:
    curDom = dom(j:j+1);
    
    % Generate a few random points to use as test values:
    x = diff(curDom) * rand(100, 1) + curDom(1);
    
    fval = feval(g, x);
    vals_exact = feval(opi{j}, x);
    err = fval - vals_exact;
    result(j) = ( norm(err-mean(err), inf) < ...
        1e7*get(f,'epslevel')*norm(vals_exact, inf) );
end
pass(11) = all( result );

%% piecewise smooth chebfun: SMOOTHFUN + SINGFUN.

% define the domain:
dom = [-1 1];
domCheck = [dom(1)+0.1 dom(2)-0.1];

op = @(x) sin(100*x)./((x-dom(1)).^0.5.*(x-dom(2)).^0.5);
f = chebfun(op, dom, 'exps', [-0.5 -0.5]);
% We temporarily disable this warning: 1
warning('off', 'CHEBFUN:SINGFUN:plus:exponentDiff');
g = cumsum(f);
warning('on', 'CHEBFUN:SINGFUN:plus:exponentDiff');

%%
% check values:

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

% g = restrict(g, domCheck);
gval = feval(g, x);
f_check = chebfun(op, domCheck, 'splitting', 'on');
g_check = cumsum(f_check);

vals_check = feval(g_check, x);
err = gval - vals_check;
err = norm(err-mean(err), inf);
tol = 5e5*get(f,'epslevel')*norm(vals_check, inf);
pass(12) = err < tol;

%% Tests for functions defined on unbounded domain:

%% Function on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];
domCheck = [-1e2 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) exp(-x.^2);
f = chebfun(op, dom);
g = cumsum(f);

gVals = feval(g, x);
opg = @(x) sqrt(pi)*erf(x)/2 + sqrt(pi)/2;
gExact = opg(x);
errg = norm(gVals - gExact, inf);
tol = 5e4*get(g,'epslevel').*get(g,'vscale');
pass(13) = errg < tol;

%% Function on [a inf]:

% Set the domain:
dom = [1 Inf];
domCheck = [1 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

% Blow-up function:
op = @(x) 5*x;
f = chebfun(op, dom, 'exps', [0 1]);
g = cumsum(f);
gVals = feval(g, x);

opg = @(x) 5*x.^2/2 - 5/2 + get(g, 'lval');
gExact = opg(x);
err = norm(gVals - gExact, inf);
tol = 100*get(g,'epslevel').*get(g,'vscale');
pass(14) = err < tol;

%% Piecewise function on [-inf b]:

% Set the domain:
dom = [-Inf -1 3*pi];
domCheck = [-1e6 -1 3*pi];

% Generate a few random points to use as test values:
x1 = diff(domCheck(1:2)) * rand(100, 1) + domCheck(1);
x2 = diff(domCheck(2:3)) * rand(100, 1) + domCheck(2);

op1 = @(x) exp(x);
op2 = @(x) sin(3*x);
f = chebfun({op1 op2}, dom);
g = cumsum(f);
g1Vals = feval(g, x1);
g2Vals = feval(g, x2);

opg1 = @(x) exp(x);
opg2 = @(x) -cos(3*x)/3 + cos(-3)/3 + exp(dom(2));
g1Exact = opg1(x1);
g2Exact = opg2(x2);
err1 = g1Vals - g1Exact;
err2 = g2Vals - g2Exact;
pass(15) = norm([err1 ; err2], inf) < 5e4*get(g,'epslevel').*get(g,'vscale');

% Check delta functions:
x = chebfun('x');
f = heaviside(x+.5)+heaviside(x) + heaviside(x-.5);
fp = diff(f);
fpp = diff(fp);
fp = cumsum(fpp);
pass(16) = norm(f - cumsum(fp)) < 100*eps;

% Bug reported by LNT #1289:
f = sign(x)-x;
f2 = cumsum(diff(f));
pass(17) = norm(f-f2) < 100*eps;

% Check quasi matrix with delta functions:
f = heaviside(x);
g = heaviside(x-.5);
s = [diff(f), diff(g)];
S = cumsum(s);
pass(18) = norm(S(:,1) - f) < 100*eps;
pass(19) = norm(S(:,2) - g) < 100*eps;


% [TODO]:  Check fractional antiderivatives once implemented.

end

function y = test_If4(x)
    y = [zeros(size(x)) zeros(size(x))];
    int1 = (-1 < x) & (x < -0.5);
    int2 = (-0.5 < x) & (x < 1);
    y(int1, 1) = sin(x(int1));
    y(int1, 2) = cos(x(int1));
    y(int2, 1) = sin(x(int2));
    y(int2, 2) = cos(x(int2));

    y = y - repmat([sin(-1) cos(-1)], size(x, 1), 1);
end

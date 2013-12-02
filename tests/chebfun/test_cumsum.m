% Test file for @chebfun/cumsum.m.

function pass = test_cumsum(pref)

% Obtain preferences.
if ( nargin == 0 )
    pref = chebpref();
end

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

% Check empty casee.
pass(1) = isempty(cumsum(chebfun()));

% Check an example with breakpoints but no impulses.
f1 = chebfun(@cos, [-1 -0.5 0.5 1], pref);
If1 = cumsum(f1);
If1_exact = @(x) sin(x) - sin(-1);
pass(2) = norm(feval(If1, xr) - If1_exact(xr), inf) < ...
    10*vscale(If1)*epslevel(If1);

% Check behavior for row chebfuns.
f1t = f1.';
If1t = cumsum(f1t);
If1t_exact = @(x) (sin(x) - sin(-1)).';
pass(3) = norm(feval(If1t, xr) - If1t_exact(xr), inf) < ...
    10*vscale(If1t)*epslevel(If1t);

% Check behavior with impulses.
f2 = f1;
f2.impulses(2, 1, 2) = 1;
f2.impulses(2, 1, 3) = 2;
If2 = cumsum(f2);
If2_exact = @test_If2;
pass(4) = norm(feval(If2, xr) - If2_exact(xr), inf) < ...
    10*vscale(If2)*epslevel(If2) && ...
    isequal(If2.impulses(:,:,2:end), [0 ; 2 ; 0 ; 0]);

% Check behavior for array-valued chebfuns.
f3 = chebfun(@(x) [cos(x) -sin(x) exp(x)], [-1 -0.5 0.5 1], pref);
If3 = cumsum(f3);
If3_exact = @(x) [sin(x) cos(x) exp(x)] - ...
    repmat([sin(-1) cos(-1) exp(-1)], size(x, 1), 1);
pass(5) = max(max(abs(feval(If3, xr) - If3_exact(xr)))) < ...
    10*vscale(If3)*epslevel(If3);

f3t = f3.';
If3t = cumsum(f3t);
If3t_exact = @(x) ([sin(x) cos(x) exp(x)] - ...
    repmat([sin(-1) cos(-1) exp(-1)], size(x, 1), 1)).';
pass(6) = max(max(abs(feval(If3t, xr) - If3t_exact(xr)))) < ...
    10*vscale(If3t)*epslevel(If3t);

f4 = chebfun(@(x) [cos(x) -sin(x)], [-1 -0.5 0.5 1], pref);
f4.impulses(2, 2, 2) = 1;
f4.impulses(1, 1, 3) = 2;
If4 = cumsum(f4);
If4_exact = @test_If4;
pass(7) = max(max(abs(feval(If4, xr) - If4_exact(xr)))) < ...
    10*vscale(If4)*epslevel(If4) && ...
    isequal(If4.impulses(:,:,2:end), [2 0 ; 0 0 ; 0 0 ; 0 0]);

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

%% Integration with singfun: 

dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

pow = -1.64;
op = @(x) (x-dom(1)).^pow;
pref.singPrefs.exponents = [pow 0];
f = chebfun(op, dom, pref);
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
g = cumsum(f);

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
    result(j) = ( norm(err-mean(err), inf) < 1e7*get(f,'epslevel')*norm(vals_exact, inf) );
end
pass(11) = all( result );

%% piecewise smooth chebfun: smoothfun + singfun & splitting on.

% define the domain:
dom = [-1 1];
domCheck = [dom(1)+0.1 dom(2)-0.1];

op = @(x) sin(50*x)./((x-dom(1)).^0.5.*(x-dom(2)).^0.5);
f = chebfun(op, dom, 'exps', [-0.5 -0.5], 'splitting', 'on');
g = cumsum(f);

% check values:

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

g = restrict(g, domCheck);
fval = feval(g, x);
f_check = chebfun(op, domCheck, 'splitting', 'on');
g_check = cumsum(f_check);

vals_check = feval(g_check, x);
err = fval - vals_check;
pass(12) = all( norm(err-mean(err), inf) < 1e3*get(f,'epslevel')*norm(vals_check, inf) );

% [TODO]:  Check fractional antiderivatives once implemented.

end

function y = test_If2(x)
    y = zeros(size(x));
    int1 = (-1 < x) & (x < -0.5);
    int2 = (-0.5 < x) & (x < 1);
    y(int1) = sin(x(int1));
    y(int2) = sin(x(int2)) + 1;

    y = y - sin(-1);
end

function y = test_If4(x)
    y = [zeros(size(x)) zeros(size(x))];
    int1 = (-1 < x) & (x < -0.5);
    int2 = (-0.5 < x) & (x < 1);
    y(int1, 1) = sin(x(int1));
    y(int1, 2) = cos(x(int1));
    y(int2, 1) = sin(x(int2));
    y(int2, 2) = cos(x(int2)) + 1;

    y = y - repmat([sin(-1) cos(-1)], size(x, 1), 1);
end

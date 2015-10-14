% Test file for 'splitting on' functionality with functions involving abs().

function pass = test_splitting_abs(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

seedRNG(6178);
pref.splitting = true;

% Some basic test functions:
FF = {@abs, @(x) abs(x).^5, @(x) abs(sin(10*x)), @(x) abs(sin(30*x))};

for j = 1:numel(FF);
    % Initialise k:
    k = 0;

    % Pick the test function:
    F = FF{j};

    % Test on [-1 1]:
    f = chebfun(F, [-1, 1], pref);
    xx = linspace(-1, 1);
    
    err = norm(feval(f, xx) - feval(F, xx), inf);
    pass(j, k+1) = err < 50*eps;
    pass(j, k+2) = err < 1000*pref.eps;

end

%% Test on singular function:

% define the domain:
dom = [-1 1];
domCheck = [dom(1)+0.1 dom(2)-0.1];

pow1 = -0.5;
pow2 = -0.2;
op = @(x) abs(cos(25*x).*((x-dom(1)).^pow1).*((dom(2)-x).^pow2));
f = chebfun(op, dom, 'exps', [pow1 pow2], 'splitting', 'on');

% check values:

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

vals_f = feval(f, x);
vals_check = feval(op, x);
err = vals_f - vals_check;

pass(j+1,:) = ( norm(err-mean(err), inf) < ...
    1e6*eps*norm(vals_check, inf) );


%% Tests for function defined on unbounded domain:

% Functions on [a inf]:
dom = [0 Inf];
domCheck = [0 100];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) 0.75+sin(10*x)./exp(x);
opAbs = @(x) abs(0.75+sin(10*x)./exp(x));
f = chebfun(op, dom, 'splitting', 'on');
g = abs(f);
gVals = feval(g, x);
gExact = opAbs(x);
err = gVals - gExact;
pass(j+2,:) = norm(err, inf) < 1e3*eps*vscale(g);


end

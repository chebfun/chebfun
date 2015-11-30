% Test file for chebfun constructor (periodic).

function pass = test_constructor_basic_periodic(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

%% Test construction of some basic functions:

% Some basic test functions:
FF = {@(x) exp(sin(pi*x)), @(x) exp([sin(pi*x), cos(pi*x)]), @(x) exp([sin(pi*x), cos(pi*x), -cos(pi*x).^2])};

for j = 1:numel(FF);
    % Initialise k:
    k = 0;

    % Pick the test function:
    F = FF{j};

    % Test on [-1 1]:
    f = chebfun(F, [-1, 1], pref, 'periodic');
    g = chebfun(F, [-1, 1], pref, 'trig');
    xx = linspace(-1, 1);
    err = norm(feval(f, xx) - F(xx), inf);
    pass(j, k+1) = (err < 10*eps*vscale(f)) && (norm(f-g) < pref.eps);
    pass(j, k+2) = err < 50*pref.eps;
    k = k + 2;

    % Test on [-1 1] (no domain passed):
    f = chebfun(F, pref, 'periodic');
    xx = linspace(-1, 1);
    err = norm(feval(f, xx) - F(xx), inf);
    pass(j, k+1) = err < 10*eps*vscale(f);
    pass(j, k+2) = err < 500*pref.eps;
    k = k + 2;

    % Test on [0 10000]:
    f = chebfun(F, [-100, 100], pref, 'periodic');
    g = chebfun(F, [-100, 100], pref, 'trig');
    xx = linspace(-100, 100);
    err = norm(feval(f, xx) - F(xx), inf);
    pass(j, k+1) = err < 1e3*eps*vscale(f);
    pass(j, k+2) = (err < 100*hscale(f)*pref.eps) && (norm(f-g) < 100*pref.eps);
    k = k + 2;
    
end


end

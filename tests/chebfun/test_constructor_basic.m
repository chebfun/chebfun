function pass = test_constructor_basic(pref)

if ( nargin == 0 )
    pref = chebfun.pref();
end

% Some basic test functions:
FF = {@sin, @(x) [sin(x), cos(x)], @(x) [sin(x), cos(x), exp(-x)]};

for j = 1:numel(FF);
    % Initialise k:
    k = 0;

    % Pick the test function:
    F = FF{j};

    % Test on [-1 1]:
    % TODO:  Remove call to cell2mat() once we've figured out the vscale issue.
    f = chebfun(F, [-1, 1], pref);
    xx = linspace(-1, 1);
    err = norm(feval(f, xx) - F(xx), inf);
    pass(j, k+1) = err < 10*epslevel(f)*vscale(f);
    pass(j, k+2) = err < 50*pref.chebfun.eps;
    k = k + 2;

    % Test on [-1 1] (no domain passed):
    % TODO:  Remove call to cell2mat() once we've figured out the vscale issue.
    f = chebfun(F, pref);
    xx = linspace(-1, 1);
    err = norm(feval(f, xx) - F(xx), inf);
    pass(j, k+1) = err < 10*epslevel(f)*vscale(f);
    pass(j, k+2) = err < 500*pref.chebfun.eps;
    k = k + 2;

    % Test on [0 10000]:
    % TODO:  Remove call to cell2mat() once we've figured out the vscale issue.
    f = chebfun(F, [0, 10000], pref);
    xx = linspace(0, 10000);
    err = norm(feval(f, xx) - F(xx), inf);
    pass(j, k+1) = err < 100*epslevel(f)*vscale(f);
    pass(j, k+2) = err < 100*hscale(f)*pref.chebfun.eps;
    k = k + 2;

    % Test on piecewise domain:
    % TODO:  Remove call to cell2mat() once we've figured out the vscale issue.
    f = chebfun(F, [-1, 0, .5, sqrt(pi/4), 1], pref);
    xx = linspace(-1, 1);
    err = norm(feval(f, xx) - F(xx), inf);
    pass(j, k+1) = err < 10*epslevel(f)*vscale(f);
    pass(j, k+2) = err < 100*pref.chebfun.eps;
    k = k + 2;
end

end

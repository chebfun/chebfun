function pass = test_splitting_abs(pref)

if ( nargin == 0 )
    pref = chebfun.pref();
end

pref.chebfun.splitting = true;

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
    pass(j, k+1) = err < 50*max(f.epslevel);
    pass(j, k+2) = err < 1000*pref.chebfun.eps;
    k = k + 2;

end

end

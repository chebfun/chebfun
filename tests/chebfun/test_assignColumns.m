% Test file for @chebfun/assignColumns.m.

function pass = test_assignColumns(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

funlist = {@chebfun, @quasimatrix};

for k = 1:2
    myfun = funlist{k};
    
    % Generate a few random points to use as test values.
    seedRNG(6178);
    x = 2 * rand(10, 1) - 1;

    % Check a few examples.
    f = myfun(@(x) [sin(x) cos(x) exp(x)], [-1 0 1], pref);
    g = myfun(@(x) x, [-1 1], pref);

    h = assignColumns(f, 1, g);
    h_exact = @(x) [x cos(x) exp(x)];
    err = feval(h, x) - h_exact(x);
    pass(k,1) = norm(err(:), inf) < 10*vscale(h)*eps;

    h = assignColumns(f, 3, g);
    h_exact = @(x) [sin(x) cos(x) x];
    err = feval(h, x) - h_exact(x);
    pass(k,2) = norm(err(:), inf) < 10*vscale(h)*eps;

    g = myfun(@(x) [x x.^2], [-1 0 1], pref);

    h = assignColumns(f, [3 1], g);
    h_exact = @(x) [x.^2 cos(x) x];
    err = feval(h, x) - h_exact(x);
    pass(k,3) = norm(err(:), inf) < 10*vscale(h)*eps;

    h = assignColumns(f, [2 2], g);
    h_exact = @(x) [sin(x) x.^2 exp(x)];
    err = feval(h, x) - h_exact(x);
    pass(k,4) = norm(err(:), inf) < 10*vscale(h)*eps;

    g = myfun(@(x) [x x.^2 x.^3], [-1 -0.5 0.5 1], pref);

    h = assignColumns(f, [1 2 3], g);
    h_exact = @(x) [x x.^2 x.^3];
    err = feval(h, x) - h_exact(x);
    pass(k,5) = norm(err(:), inf) < 10*vscale(h)*eps;

    hc = assignColumns(f, ':', g);
    pass(k,6) = isequal(hc, h);

    h  = assignColumns(f, [2 1], [-0.5 0.5]);
    h_exact = @(x) [(0.5 + 0*x) (-0.5 + 0*x) exp(x)];
    err = feval(h, x) - h_exact(x);
    pass(k,7) = norm(err(:), inf) < 10*vscale(h)*eps;

    h = assignColumns(f.', [2 1], [-0.5 ; 0.5]);
    h_exact = @(x) [(0.5 + 0*x) (-0.5 + 0*x) exp(x)].';
    err = feval(h, x) - h_exact(x);
    pass(k,8) = norm(err(:), inf) < 10*vscale(h)*eps;

    % Check error conditions.
    try
        h = assignColumns(f, [1 2 3], g.');
        pass(k,9) = false;
    catch ME
        pass(k,9) = strcmp(ME.identifier, ...
            'CHEBFUN:CHEBFUN:assignColumns:numCols');
    end

    try
        h = assignColumns(f, [1 2], g);
        pass(k,10) = false;
    catch ME
        pass(k,10) = strcmp(ME.identifier, ...
            'CHEBFUN:CHEBFUN:assignColumns:numCols');
    end

    try
        g2 = myfun(@(x) [x x.^2 x.^3], [0 1], pref);
        h = assignColumns(f, [1 2 3], g2);
        pass(k,11) = false;
    catch ME
        pass(k,11) = strcmp(ME.identifier, ...
            'CHEBFUN:CHEBFUN:assignColumns:domain');
    end

    % Test assign outside of dimension of f:
    h = assignColumns(f, [1 2 4], g);
    h_exact = @(x) [x x.^2 exp(x) x.^3];
    err = feval(h, x) - h_exact(x);
    pass(k,12) = norm(err(:), inf) < 10*vscale(h)*eps;

    %% Test for function defined on unbounded domain:
    
    % Functions on [-inf b]:
    
    % Set the domain:
    dom = [-Inf -3*pi];
    domCheck = [-1e6 -3*pi];
    
    % Generate a few random points to use as test values:
    x = diff(domCheck) * rand(100, 1) + domCheck(1);
    
    % Array-valued function:
    op = @(x) [exp(x) x.*exp(x) (1-exp(x))./x];
    opg = @(x) exp(-x.^2);
    oph = @(x) [exp(x) exp(-x.^2) (1-exp(x))./x];
    
    f = myfun(op, dom);
    g = myfun(opg, dom);
    h = assignColumns(f, 2, g);
    err = feval(h, x) - oph(x);    
    pass(k,13) = norm(err, inf) < 1e2*vscale(h)*eps;
       
 
    %% Test removing columns (#889)
    
    x = chebfun('x');
    xx = [x 2*x];
    if ( k == 2 )
        xx = cheb2quasi(xx);
    end
    xx(:,2) = [];
    pass(k,14) = numColumns(xx) == 1 && norm(x - xx) < eps;
    
    x = chebfun('x');
    xx = [x 2*abs(x)];
    if ( k == 1 )
        xx = quasi2cheb(xx);
    end
    xx(:,2) = [];
    pass(k,15) = numColumns(xx) == 1 && norm(x - xx) < eps;
    
end

end

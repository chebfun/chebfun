% Test file for @chebfun/assignColumns.m.

function pass = test_assignColumns(pref)

if ( nargin < 1 )
    pref = chebpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

funlist = {@chebfun, @quasimatrix};

for k = 1:2
    myfun = funlist{k};


    % Check a few examples.
    f = myfun(@(x) [sin(x) cos(x) exp(x)], [-1 0 1], pref);
    g = myfun(@(x) x, [-1 1], pref);

    h = assignColumns(f, 1, g);
    h_exact = @(x) [x cos(x) exp(x)];
    err = feval(h, x) - h_exact(x);
    pass(k,1) = norm(err(:), inf) < 10*vscale(h)*epslevel(h);

    h = assignColumns(f, 3, g);
    h_exact = @(x) [sin(x) cos(x) x];
    err = feval(h, x) - h_exact(x);
    pass(k,2) = norm(err(:), inf) < 10*vscale(h)*epslevel(h);

    g = myfun(@(x) [x x.^2], [-1 0 1], pref);

    h = assignColumns(f, [3 1], g);
    h_exact = @(x) [x.^2 cos(x) x];
    err = feval(h, x) - h_exact(x);
    pass(k,3) = norm(err(:), inf) < 10*vscale(h)*epslevel(h);

    h = assignColumns(f, [2 2], g);
    h_exact = @(x) [sin(x) x.^2 exp(x)];
    err = feval(h, x) - h_exact(x);
    pass(k,4) = norm(err(:), inf) < 10*vscale(h)*epslevel(h);

    g = myfun(@(x) [x x.^2 x.^3], [-1 -0.5 0.5 1], pref);

    h = assignColumns(f, [1 2 3], g);
    h_exact = @(x) [x x.^2 x.^3];
    err = feval(h, x) - h_exact(x);
    pass(k,5) = norm(err(:), inf) < 10*vscale(h)*epslevel(h);

    hc = assignColumns(f, ':', g);
    pass(k,6) = isequal(hc, h);

    h  = assignColumns(f, [2 1], [-0.5 0.5]);
    h_exact = @(x) [(0.5 + 0*x) (-0.5 + 0*x) exp(x)];
    err = feval(h, x) - h_exact(x);
    pass(k,7) = norm(err(:), inf) < 10*vscale(h)*epslevel(h);

    h = assignColumns(f.', [2 1], [-0.5 ; 0.5]);
    h_exact = @(x) [(0.5 + 0*x) (-0.5 + 0*x) exp(x)].';
    err = feval(h, x) - h_exact(x);
    pass(k,8) = norm(err(:), inf) < 10*vscale(h)*epslevel(h);

    % Check error conditions.
    try
        h = assignColumns(f, [1 2 3], g.');
        pass(k,9) = false;
    catch ME
        pass(k,9) = strcmp(ME.identifier, 'CHEBFUN:assignColumns:numCols');
    end

    try
        h = assignColumns(f, [1 2], g);
        pass(k,10) = false;
    catch ME
        pass(k,10) = strcmp(ME.identifier, 'CHEBFUN:assignColumns:numCols');
    end

    try
        g2 = myfun(@(x) [x x.^2 x.^3], [0 1], pref);
        h = assignColumns(f, [1 2 3], g2);
        pass(k,11) = false;
    catch ME
        pass(k,11) = strcmp(ME.identifier, 'CHEBFUN:assignColumns:domain');
    end

    try
        h = assignColumns(f, [1 2 4], g);
        pass(k,12) = false;
    catch ME
        pass(k,12) = strcmp(ME.identifier, 'CHEBFUN:assignColumns:dims');
    end
    
end

end

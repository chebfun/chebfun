function pass = test_extractBoundaryRoots(pref)

if ( nargin < 1 )
    pref = chebtech.techPref();
end

% Generate random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else
        testclass = chebtech2();
    end
    
    %% Test roots at left endpoint:
    ml = 3;
    f = testclass.make(@(x) sin(2*x).*((1 + x).^ml), [], pref);
    [g, l, r] = extractBoundaryRoots(f);
    gexact = testclass.make(@(x) sin(2*x), [], pref);
    err = feval(g, x) - feval(gexact, x);
    pass(n, 1) = (norm(err, Inf) < (1e2^ml)*f.epslevel) && (l == ml);
    
    %% Test roots at right endpoint:
    mr = 2;
    f = testclass.make(@(x) sin(cos(3*x)).*((1 - x).^mr), [], pref);
    [g, l, r] = extractBoundaryRoots(f);
    gexact = testclass.make(@(x) sin(cos(3*x)), [], pref);
    err = feval(g, x) - feval(gexact, x);
    pass(n, 2) = (norm(err, Inf) < (5e2^mr)*f.epslevel) && (r == mr);
    
    %% Test roots at both endpoints:
    ml = 1;
    mr = 2;
    f = testclass.make(@(x) exp(x).*((1 + x).^ml).*((1 - x).^mr), [], pref);
    [g, l, r] = extractBoundaryRoots(f);
    gexact = testclass.make(@(x) exp(x), [], pref);
    err = feval(g, x) - feval(gexact, x);
    pass(n, 3) = (norm(err, Inf) < (1e1^(ml + mr))*f.epslevel) && ...
        (l == ml) && (r == mr);
    
    %% Test complex case:
    ml = 1;
    mr = 2;
    f_op = @(x) (x.^2 + exp(x) + 1i*cos(2*x)).*((1 + x).^ml).*((1 - x).^mr);
    f = testclass.make(f_op, [], pref);
    [g, l, r] = extractBoundaryRoots(f);
    gexact = testclass.make(@(x) x.^2 + exp(x) + 1i*cos(2*x), [], pref);
    err = feval(g, x) - feval(gexact, x);
    pass(n, 4) = (norm(err, Inf) < (2e1^(ml + mr))*f.epslevel) && ...
        (l == ml) && (r == mr);
    
    %% Test when no roots:
    f = testclass.make(@(x) sin(1 - x)./(1 - x), [], pref);
    [g, l, r] = extractBoundaryRoots(f);
    err = feval(g, x) - feval(f, x);
    pass(n, 5) = (norm(err, Inf) == 0) && (l == 0) && (r == 0);
    
    %% Test when roots are not explicit:
    f = testclass.make(@(x) sin(1 - x), [], pref);
    [g, l, r] = extractBoundaryRoots(f);
    gexact = testclass.make(@(x) sin(1 - x)./(1 - x), [], pref);
    err = feval(g, x) - feval(gexact, x);
    pass(n, 6) = (norm(err, Inf) < 1e2*f.epslevel) && (r == 1);
    
    %% Test array-valued case:
    f_op = @(x) [sin(x).*((1-x).^2) cos(x.^2).*(1+x).*(1-x)];
    f = testclass.make(f_op, [], pref);
    [g, l, r] = extractBoundaryRoots(f);
    gexact = testclass.make(@(x) [sin(x) cos(x.^2)], [], pref);
    err = feval(g, x) - feval(gexact, x);
    pass(n, 7) = (norm(err, Inf) < (1e2^2)*max(f.epslevel)) && all(l == [0 1]) && ...
        all(r == [2 1]);
    
    %% Test on full arguments:
    ml = 1;
    mr = 2;
    f = testclass.make(@(x) exp(x).*((1 + x).^ml).*((1 - x).^mr), [], pref);
    g = extractBoundaryRoots(f, [ml; mr]);
    gexact = testclass.make(@(x) exp(x), [], pref);
    err = feval(g, x) - feval(gexact, x);
    pass(n, 8) = (norm(err, Inf) < (1e1^(ml + mr))*f.epslevel);
    
    %% Test on wrong multiplicities supplied by users:
    ml = 1;
    mr = 2;
    f = testclass.make(@(x) exp(x).*((1 + x).^ml).*((1 - x).^mr), [], pref);
    [g, l, r] = extractBoundaryRoots(f, [ml+1; mr+2]);
    gexact = testclass.make(@(x) exp(x), [], pref);
    err = feval(g, x) - feval(gexact, x);
    pass(n, 9) = ( (norm(err, Inf) < (1e1^(ml + mr))*f.epslevel) && ...
        ( l == ml ) && ( r == mr ) );    
    
end

end

function pass = test_sign(pref)

if ( nargin == 0 ) 
    pref = chebtech.techPref();
end

for type = 1:2
    
    if ( type == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end
    
    % Test a positive function:
    F = @(x) sin(x) + 2;
    f = testclass.make(@(x) F(x), [], pref);
    h = sign(f);
    pass(type, 1) = normest(h - 1) < 10*f.epslevel;
    
    % Test a negative function:
    f2 = testclass.make(@(x) -F(x), [], pref);
    h = sign(f2);
    pass(type, 2) = normest(h + 1) < 10*f.epslevel;
    
    % Test a complex-valued function:
    F = @(x) exp(1i*pi*x);
    f = testclass.make(@(x) F(x), [], pref);
    h = sign(f);
    pass(type, 3) = normest(h - f) < 10*f.epslevel;
    
    % Test a complex array-valued function:
    xx = linspace(-.95, .97);
    F = @(x) [(2+sin(x)).*exp(1i*pi*x), -(2+sin(x)).*exp(1i*pi*x), 2+sin(x)];
    f = testclass.make(@(x) F(x), [], pref);
    ff = feval(f, xx);
    gg = ff./abs(ff);
    h = sign(f);
    hh = feval(h, xx);
    pass(type, 4) = norm(hh - gg, inf) < 10*max(f.epslevel);
    
end

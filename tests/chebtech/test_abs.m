function pass = test_abs(pref)

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
    h = abs(f);
    pass(type, 1) = normest(h - f) < 10*eps;
    
    % Test a negative function:
    f2 = testclass.make(@(x) -F(x), [], pref);
    h = abs(f2);
    pass(type, 2) = normest(h + f2) < 10*eps;
    
    % Test a complex-valued function:
    F = @(x) exp(1i*pi*x);
    f = testclass.make(@(x) F(x), [], pref);
    h = abs(f);
    pass(type, 3) = normest(h - 1) < 1e2*eps;
    
    % Test a complex array-valued function:
    F = @(x) [(2+sin(x)).*exp(1i*pi*x), -(2+sin(x)).*exp(1i*pi*x), 2+sin(x)];
    f = testclass.make(@(x) F(x), [], pref);
    g = testclass.make(@(x) [2+sin(x), 2+sin(x), 2+sin(x)]);
    h = abs(f);
    pass(type, 4) = normest(h - g) < 1e1*eps;
    
end

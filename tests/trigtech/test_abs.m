function pass = test_abs(pref)

if ( nargin == 0 ) 
    pref = trigtech.techPref();
end

testclass = trigtech();
    
% Test a positive function:
F = @(x) sin(pi*x) + 2;
f = testclass.make(@(x) F(x), [], pref);
h = abs(f);
pass(1) = normest(h - f) < 10*f.epslevel;

% Test a negative function:
f2 = testclass.make(@(x) -F(x), [], pref);
h = abs(f2);
pass(2) = normest(h + f2) < 10*f.epslevel;

% Test a complex-valued function:
F = @(x) exp(1i*pi*x);
f = testclass.make(@(x) F(x), [], pref);
h = abs(f);
pass(3) = normest(h - 1) < 10*f.epslevel;

% Test a complex array-valued function:
F = @(x) [(2+sin(pi*x)).*exp(1i*pi*x), -(2+sin(pi*x)).*exp(1i*pi*x), 2+sin(pi*x)];
f = testclass.make(@(x) F(x), [], pref);
g = testclass.make(@(x) [2+sin(pi*x), 2+sin(pi*x), 2+sin(pi*x)]);
h = abs(f);
pass(4) = normest(h - g) < 10*max(f.epslevel);

end

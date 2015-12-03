% Test file for trigtech/sign.m

function pass = test_sign(pref)

if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();
    
% Test a positive function:
F = @(x) sin(pi*x) + 2;
f = testclass.make(@(x) F(x), [], pref);
h = sign(f);
pass(1) = normest(h - 1) < 10*eps;

% Test a negative function:
f2 = testclass.make(@(x) -F(x), [], pref);
h = sign(f2);
pass(2) = normest(h + 1) < 10*eps;

% Test a complex-valued function:
F = @(x) exp(1i*pi*x);
f = testclass.make(@(x) F(x), [], pref);
h = sign(f);
pass(3) = normest(h - f) < 10*eps;

% Test a complex array-valued function:
xx = linspace(-.95, .97);
F = @(x) [(2+sin(pi*x)).*exp(1i*pi*x), -(2+sin(pi*x)).*exp(1i*pi*x), 2+sin(pi*x)];
f = testclass.make(@(x) F(x), [], pref);
ff = feval(f, xx);
gg = ff./abs(ff);
h = sign(f);
hh = feval(h, xx);
pass(4) = norm(hh - gg, inf) < 10*eps;
    
end

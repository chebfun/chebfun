% Test file for fourtech/length.m

function pass = test_length(pref)

if ( nargin < 1 )
    pref = fourtech.techPref();
end

testclass = fourtech();

f = testclass.make(@(x) tanh(sin(pi*x)), [], pref);
pass(1) = length(f) == size(f.coeffs, 1);

f = testclass.make(@(x) tanh([sin(pi*x), cos(pi*x), 1i*exp(x)]), [], pref);
pass(2) = length(f) == size(f.coeffs, 1);

p = pref;
p.fixedLength = 101;
f = testclass.make(@(x) [sin(pi*x), cos(pi*x), exp(1i*pi*x)], [], p);
pass(3) = length(f) == 101;

end

% Test file for trigtech/quadwts.m

function pass = test_quadpts(pref)

if ( nargin == 1 )
    pref = trigtech.techPref();
end

testTech = trigtech();
n = 10;
w = testTech.quadwts(n);

% Test against some low order trig functions:
x = testTech.trigpts(n);
pass(1) = abs(sum(w) - 2) < 2*eps;
pass(2) = abs(w*sin(pi*x)) < eps;
pass(3) = abs(w*(sin(2*pi*x).*cos(2*pi*x))) < 2*eps;
pass(4) = abs(w*(sin(2*pi*x).^2)-1) < 2*eps;
pass(5) = abs(w*(cos(2*pi*x).^2)-1) < 2*eps;

% Test against some basic properties:
pass(6) = isempty(testTech.quadwts(0));
pass(7) = testTech.quadwts(1) == 2;
pass(8) = all(testTech.quadwts(2) == 1);
    
end

    

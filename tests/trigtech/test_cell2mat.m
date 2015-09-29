% Test file for trigtech/cell2mat.m

function pass = test_cell2mat(pref)

if ( nargin < 2 )
    pref = trigtech.techPref();
end

testclass = trigtech();

f = testclass.make(@(x) [sin(pi*x) cos(pi*x) exp(2*1i*pi*x)], [], pref);
g = testclass.make(@(x) sin(pi*x), [], pref);
h = testclass.make(@(x) [cos(pi*x) exp(2*1i*pi*x)], [], pref);

F = cell2mat([g h]);
pass(1) = all( sum(F - f) < max(vscale(f)*eps) );

end

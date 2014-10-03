% Test file for trigtech/mat2cell.m

function pass = test_mat2cell(pref)

if ( nargin < 2 )
    pref = trigtech.techPref();
end

testclass = trigtech();

f = testclass.make(@(x) [sin(pi*x) cos(pi*x) exp(cos(pi*x))], [], pref);
g = testclass.make(@(x) sin(pi*x), [], pref);
h = testclass.make(@(x) [cos(pi*x) exp(cos(pi*x))], [], pref);

F = mat2cell(f, 1, [1 2]);
pass(1) = sum(F{1} - g) < 10*g.vscale.*g.epslevel;
pass(2) = all( sum(F{2} - h) < 10*max(h.vscale.*h.epslevel) );

end

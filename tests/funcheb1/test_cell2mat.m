function pass = test_cell2mat(pref)

if ( nargin < 2 )
    pref = funcheb.pref;
end

f = funcheb1(@(x) [sin(x) cos(x) exp(x)], pref);
g = funcheb1(@(x) sin(x), pref);
h = funcheb1(@(x) [cos(x) exp(x)], pref);

F = cell2mat([g h]);
pass(1) = all( sum(F - f) < f.vscale*f.epslevel );

end

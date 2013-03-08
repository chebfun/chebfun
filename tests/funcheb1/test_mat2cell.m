function pass = test_mat2cell(pref)

if ( nargin < 2 )
    pref = funcheb1.pref;
end

f = funcheb1(@(x) [sin(x) cos(x) exp(x)], pref);
g = funcheb1(@(x) sin(x), pref);
h = funcheb1(@(x) [cos(x) exp(x)], pref);

F = mat2cell(f, 1, [1 2]);
pass(1) = sum(F(1) - g) < g.epslevel;
pass(2) = all( sum(F(1) - g) < h.vscale*h.epslevel );

end
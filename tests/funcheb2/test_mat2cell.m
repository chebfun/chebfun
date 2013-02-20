function pass = test_mat2cell(pref)

if ( nargin < 2 )
    pref = funcheb2.pref;
end

f = funcheb2(@(x) [sin(x) cos(x) exp(x)], pref);
g = funcheb2(@(x) sin(x), pref);
h = funcheb2(@(x) [cos(x) exp(x)], pref);

F = mat2cell(f, 1, [1 2]);
pass(1) = sum(F(1) - g) < g.epslevel;
pass(2) = all( sum(F(1) - g) < h.vscale*h.epslevel );
pass(3) = F(1).vscale == sin(1);
pass(4) = all( F(2).vscale == [cos(0) exp(1)] );

end
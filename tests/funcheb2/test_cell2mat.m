function pass = test_cell2mat(pref)

if ( nargin < 2 )
    pref = funcheb2.pref;
end

f = funcheb2(@(x) [sin(x) cos(x) exp(x)], pref);
g = funcheb2(@(x) sin(x), pref);
h = funcheb2(@(x) [cos(x) exp(x)], pref);

F = cell2mat([g h]);
pass(1) = all( sum(F - f) < f.vscale*f.epslevel );
pass(2) = all( F.vscale == [sin(1) cos(0) exp(1)] );

end

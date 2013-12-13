function pass = test_extractColumns(pref)

if ( nargin == 0 )
    pref = chebpref();
end

f = chebfun(@(x) [sin(x), cos(x), exp(x)], pref);
g = chebfun(@(x) [sin(x), cos(x)], pref);
h = extractColumns(f, 1:2);
pass(1) = size(h, 2) == 2 && normest(g - h) < epslevel(f) && ...
    norm(g.pointValues - h.pointValues, inf) < epslevel(f);
h = f(:,1:2);
pass(2) = size(h, 2) == 2 && normest(g - h) < epslevel(f) && ...
    norm(g.pointValues - h.pointValues, inf) < epslevel(f);

g = chebfun(@(x) [sin(x), sin(x), exp(x), cos(x)], pref);
h = extractColumns(f, [1 1 3 2]);
pass(3) = size(h, 2) == 4 && normest(g - h) < epslevel(f) && ...
    norm(g.pointValues - h.pointValues, inf) < epslevel(f);

end


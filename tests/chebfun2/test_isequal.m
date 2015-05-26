% Test file for @chebfun2/isequal.m

function pass = test_isequal(pref)

if ( nargin < 1 )
    pref = chebfunpref;
end

f = chebfun2(@(x, y) cos(x.*y));
g = chebfun2(@(x, y) sin(x + y.^2));

pass(1) = isequal(f, f);
pass(2) = ~isequal(f, g);

end

function pass = test_chebpolyvalm(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

seedRNG(0)
p_cheb = rand(3, 1) + rand(3, 1)*1i;
A = rand(3) + rand(3)*1i;
f = chebtech2({[], p_cheb});
p_poly = poly(f);

err = norm(polyvalm(p_poly, A) - chebpolyvalm(flipud(p_cheb), A));
tol = 100*pref.chebfuneps;
pass(1) = err < tol;

err = norm(polyvalm(p_poly.', A) - chebpolyvalm((flipud(p_cheb)).', A));
tol = 100*pref.chebfuneps;
pass(2) = err < tol;

try
    chebpolyvalm(p_cheb, rand(3, 4))
    pass(3) = false;
catch ME
    pass(3) = strcmp(ME.identifier, 'CHEBFUN:chebpolyvalm:square');
end

try
    chebpolyvalm(rand(2, 3), A)
    pass(4) = false;
catch ME
    pass(4) = strcmp(ME.identifier, 'CHEBFUN:chebpolyvalm:vector');
end

end
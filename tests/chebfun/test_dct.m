function pass = test_dct(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% TODO: include more extensive tests.

n = 5;
seedRNG(42)
R = rand(n);

tol = n*eps;

pass(1) = norm(chebfun.idct(chebfun.dct(R,1),1) - R) < tol;
pass(2) = norm(chebfun.idct(chebfun.dct(R,2),2) - R) < tol;
pass(3) = norm(chebfun.idct(chebfun.dct(R,3),3) - R) < tol;
pass(4) = norm(dct(R) - chebfun.dct(R,2)) < tol;

end

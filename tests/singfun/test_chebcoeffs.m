function pass = test_chebpoly(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

f = singfun(@(x) sqrt(1-x));
c = chebcoeffs(f, 10);

exact = [-2*sqrt(2)/35, sqrt(2)];
pass(1) = length(c) == 10;
pass(2) = norm(c([7 10]) - exact) < 10*get(f, 'epslevel');

end
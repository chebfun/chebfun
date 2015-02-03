function pass = test_chebcoeffs(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

f = singfun(@(x) sqrt(1-x));
c = chebcoeffs(f, 10);

exact = [sqrt(2) ; -2*sqrt(2)/35]*(2/pi);
pass(1) = length(c) == 10;
pass(2) = norm(c([1 ; 4]) - exact) < 10*get(f, 'epslevel');

end
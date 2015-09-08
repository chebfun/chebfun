function pass = test_chebcoeffs(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

f = singfun(@(x) sqrt(1-x));
c = chebcoeffs(f, 10);

exact = [sqrt(2) ; -2*sqrt(2)/35]*(2/pi);
pass(1) = length(c) == 10;
pass(2) = norm(c([1 ; 4]) - exact) < 10*eps;

% Check second-kind coefficients.
f = singfun(@(x) (1 + x.^2 + x.^3).*sqrt(1 - x.^2), ...
    struct('exponents', [0.5 0.5]));
c = chebcoeffs(f, 10, 2);

exact = [8/5 ; 16/315]*(2/pi);
pass(3) = length(c) == 10;
pass(4) = norm(c([1 ; 4]) - exact) < 10*eps;

end

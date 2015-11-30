function pass = test_var(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

f = chebfun({@(x) exp(4*pi*1i*x), @exp, @exp}, [-1 0 0.5 1], pref);
pass(1) = abs(var(f) - exp(1)/2) < 10*vscale(f)*eps;

ft = f.';
pass(2) = abs(var(ft) - exp(1)/2) < 10*vscale(ft)*eps;

f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1], pref);
varf_exact = [(1 - sin(2)/2)/2 ((1 + sin(1)*cos(1))/2 - sin(1)^2) ...
    (1 - exp(-2))/2];
pass(3) = norm(var(f) - varf_exact, inf) < 10*vscale(f)*eps;

end

function pass = test_chebcoeffsy(pref)

if ( nargin == 0 ) 
    pref = chebfunpref();
end

% Test on piecewise-smooth chebfun
f = chebfun(@(x)sin(100*x), 'splitting', 'on');
n = 20;
p = chebcoeffs(f, n);

g = chebfun(@(x)sin(100*x));
c = g.funs{1}.onefun.coeffs;
cc = flipud(c(end-n+1:end));

err = cc-p;
pass(1) = norm(err, inf) < vscale(f).*epslevel(f);

end
function pass = test_chebpoly(pref)

if ( nargin == 0 ) 
    pref = chebpref();
end

% Test on piecewise-smooth chebfun
f = chebfun(@(x)sin(100*x), 'splitting', 'on');
n = 20;
p = chebpoly(f, nan, n).';

g = chebfun(@(x)sin(100*x));
c = g.funs{1}.onefun.coeffs;
cc = c(end-n+1:end);

err = cc-p;
pass(1) = norm(err, inf) < vscale(f)*epslevel(f);

end
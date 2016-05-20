function pass = test_get(pref)
% Test GET.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e2 * pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) cos(x.*y.*z)); 
[core, cols, rows, tubes] = tucker(f);

pass(1) = norm([-1 1 -1 1 -1 1] - f.domain) < tol;

pass(2) = norm(core(:) - f.core(:)) < tol;

pass(3) = norm(cols - f.cols) < tol;

pass(4) = norm(rows - f.rows) < tol;

pass(5) = norm(tubes - f.tubes) < tol;

end
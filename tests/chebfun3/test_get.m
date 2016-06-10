function pass = test_get(pref)
% Test GET.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e2 * pref.cheb3Prefs.chebfun3eps;
j = 1; 

f = chebfun3(@(x,y,z) cos(x.*y.*z)); 
[core, cols, rows, tubes] = tucker(f);
pass(j) = norm([-1 1 -1 1 -1 1] - f.domain) < tol; j = j + 1; 
pass(j) = norm(core(:) - f.core(:)) < tol; j = j + 1; 
pass(j) = norm(cols - f.cols) < tol; j = j + 1; 
pass(j) = norm(rows - f.rows) < tol; j = j + 1; 
pass(j) = norm(tubes - f.tubes) < tol; j = j + 1; 

end
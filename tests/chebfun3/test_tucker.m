function pass = test_tucker(pref)
% Test ST, i.e., Slice-Tucker decomposition 

if ( nargin == 0) 
    pref = chebfunpref; 
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) cos(x.*y.*z));
[fCore, fCols, fRows, fTubes] = tucker(f); 
x = linspace(-1,1)';
[xx,yy,zz] = ndgrid(x);
fVals = f(xx, yy, zz);
stVal = chebfun3.txm(chebfun3.txm(chebfun3.txm(fCore, fCols(x, :), 1), ...
    fRows(x, :), 2), fTubes(x, :), 3);
err = norm(stVal(:) - fVals(:));
pass(1) = err < tol;

f = chebfun3(@(x,y,z) cos(x.*y.*z), [-3 4 -1 3 2 4]);
[fCore, fCols, fRows, fTubes] = tucker(f); 
x = linspace(-3,4)';
y = linspace(-1,3)';
z = linspace(2,4)';
[xx,yy, zz] = ndgrid(x, y, z);
fVals = f(xx, yy, zz);
stVal = chebfun3.txm(chebfun3.txm(chebfun3.txm(fCore, fCols(x, :), 1), ...
    fRows(y, :), 2), fTubes(z, :), 3);
err = norm(stVal(:) - fVals(:));
pass(2) = err < tol;

core = tucker(f); 
pass(3) = norm(fCore(:) - core(:)) < tol;

end
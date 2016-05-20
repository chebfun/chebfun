function pass = test_minus(pref)
% Test for chebfun3/minus.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e5 * pref.cheb3Prefs.chebfun3eps;

dom = [-1 1 -1 1 -1 1; 
       -2 2 -2 2 -2 2; 
       -1 pi 0 2*pi -pi pi];
   
ff = @(x,y,z) cos(x.*y.*z);
gg = @(x,y,z) x + y + z + x.*y.*z;

for j = 1:3
    f = chebfun3(ff, dom(j,:));
    g = chebfun3(gg, dom(j,:));

    FminusG = chebfun3(@(x,y,z) ff(x,y,z) - gg(x, y, z), dom(j,:));
    
    tolk = norm(dom(j, :), inf) * tol;
    pass(j) = norm((f-g) - FminusG) < 10*tolk;
end

end
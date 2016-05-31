function pass = test_minus(pref)
% Test for chebfun3/minus.
% [reviewed by LNT 31.05.16]

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.cheb3Prefs.chebfun3eps;

dom = [-1 1 -1 1 -1 1; 
       -1 pi 0 2*pi -pi pi];
   
ff = @(x,y,z) cos(x.*y.*z);
gg = @(x,y,z) x + y + z + x.*y.*z;

for j = 1:2
    f = chebfun3(ff, dom(j,:));
    g = chebfun3(gg, dom(j,:));

    FminusG = chebfun3(@(x,y,z) ff(x,y,z) - gg(x, y, z), dom(j,:));
    
    tolj = norm(dom(j, :), inf) * tol;
    pass(j) = norm((f-g) - FminusG) < tolj;
end

end
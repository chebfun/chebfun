function pass = test_minus( pref ) 
% Test for chebfun3/minus.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e5 * pref.cheb3Prefs.chebfun3eps;
j = 1;

dom = [-1 1 -1 1 -1 1; -2 2 -2 2 -2 2; -1 pi 0 2*pi -pi pi];
ff = @(x,y,z) cos(x.*y.*z);
gg = @(x,y,z) x + y + z + x.*y.*z;

for k = 1:size(dom, 1)
    f = chebfun3(ff, dom(k,:));
    g = chebfun3(gg, dom(k,:));

    FminusG = chebfun3(@(x,y,z) ff(x,y,z) - gg(x, y, z), dom(k,:));
    
    tolk = norm(dom(k, :), inf) * tol;
    pass(j) = norm((f-g) - FminusG) < 10*tolk; 
    j = j + 1;
end

end
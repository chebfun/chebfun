function pass = test_uplus(pref)
% Test chebfun3/uplus.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e5 * pref.cheb3Prefs.chebfun3eps;

dom = [-1 1 -1 1 -1 1; 
       -2 2 -2 2 -2 2; 
       -1 pi 0 2*pi -pi pi];

for j = 1 : size(dom, 1)
    f = chebfun3(@(x,y,z) cos(x.*y.*z), dom(j,:));
    
    uplusF = f;
    tolk = norm(dom(j, :), inf) * tol;
    pass(j) = norm(f - uplusF) < tolk;
end

end
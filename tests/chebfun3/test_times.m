function pass = test_times(pref)
% Test chebfun3/times

if ( nargin == 0) 
    pref = chebfunpref; 
end

tol = 1000*pref.cheb3Prefs.chebfun3eps;
j = 1; 

f = chebfun3(@(x,y,z) cos(x.*y.*z)); 
h = chebfun3(@(x,y,z) 2*cos(x.*y.*z)); 
k = chebfun3(@(x,y,z) cos(x.*y.*z).^2); 

pass(j) = norm(f.*2 - h) < tol; j = j + 1; 
pass(j) = norm(f*2 - h) < tol; j = j + 1; 
pass(j) = norm(2*f - h) < tol; j = j + 1; 
pass(j) = norm(2.*f - h) < tol; j = j + 1; 
pass(j) = norm(f.^2 - k) < tol; j = j + 1; 
pass(j) = norm(f.*f - k) < tol; j = j + 1;

ff = @(x,y,z) cos(x.*y.*z);
gg = @(x,y,z) x + y + z + x.*y.*z;
dom = [-1 1 -1 1 -1 1; 
       -2 2 -2 2 -2 2; 
       -1 pi 0 2*pi -pi pi];
for r = 1:size(dom, 1)
    f = chebfun3(ff, dom(r,:));
    g = chebfun3(gg, dom(r,:));
    FtimesG = chebfun3(@(x,y,z) ff(x,y,z).*gg(x, y, z), dom(r,:));
    tolk = norm(dom(r, :), inf) * tol;
    pass(j) = norm( f.*g - FtimesG ) < 10*tolk; 
    j = j + 1;
end

end
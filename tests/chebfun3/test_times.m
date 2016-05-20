function pass = test_times(pref)
% Test chebfun3/times

if ( nargin == 0) 
    pref = chebfunpref; 
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) cos(x.*y.*z)); 
h = chebfun3(@(x,y,z) 2*cos(x.*y.*z)); 
k = chebfun3(@(x,y,z) cos(x.*y.*z).^2); 

pass(1) = norm(f.*2 - h) < tol;
pass(2) = norm(f*2 - h) < tol;
pass(3) = norm(2*f - h) < tol;
pass(4) = norm(2.*f - h) < tol;
pass(5) = norm(f.^2 - k) < tol;
pass(6) = norm(f.*f - k) < tol;

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
    pass(6+r) = norm(f.*g - FtimesG) < 10*tolk;
end

end
function pass = test_times(pref)
% Test chebfun3/times

if ( nargin == 0) 
    pref = chebfunpref; 
end
tol = 1e4*pref.cheb3Prefs.chebfun3eps;

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
dom = [-1 1 -1 1 -1 1];
f = chebfun3(ff, dom);
g = chebfun3(gg, dom);
FtimesG = chebfun3(@(x,y,z) ff(x,y,z).*gg(x, y, z), dom);
tolj = norm(dom, inf) * tol;
pass(7) = norm(f.*g - FtimesG) < tolj;

dom = [-2 2 -2 2 -2 2];
f = chebfun3(ff, dom);
g = chebfun3(gg, dom);
FtimesG = chebfun3(@(x,y,z) ff(x,y,z).*gg(x, y, z), dom);
tolj = norm(dom, inf) * tol;
pass(8) = norm(f.*g - FtimesG) < tolj;

dom = [0 pi 0 pi -pi/2 pi/2];
f = chebfun3(ff, dom);
g = chebfun3(gg, dom);
FtimesG = chebfun3(@(x,y,z) ff(x,y,z).*gg(x, y, z), dom);
tolj = norm(dom, inf) * tol;
pass(9) = norm(f.*g - FtimesG) < tolj;

end
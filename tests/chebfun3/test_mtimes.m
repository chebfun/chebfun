function pass = test_mtimes(pref)
% Test chebfun3/mtimes
 
if ( nargin == 0) 
    pref = chebfunpref; 
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;
fDom = [-1 1 0 2 -4 -2];
f = chebfun3(@(x,y,z) x + cos(y+2*z), fDom);
 
% Apply f to a scalar:
c = 10;
h = f*c;
hExact = chebfun3(@(x,y,z) c.*(x + cos(y+2*z)), fDom);
pass(1) = norm(hExact - h) < tol;
 
% Apply f to a 1D Chebfun:
g1D = chebfun(@(x) 2*x+3);
h = f*g1D;
hExact = chebfun2(@(y,z) 4/3+6*cos(y+2*z), fDom(3:6));
pass(2) = norm(hExact - h) < tol;
 
% Apply f to a Chebfun2:
gDom = [-1 1 3 5];
g2D = chebfun2(@(x,t) x+t, gDom);
h = f*g2D;
hExact = chebfun3(@(t,y,z) 2/3+2*t.*cos(y+2*z), [gDom(3:4) fDom(3:6)]);
pass(3) = norm(hExact - h) < tol;
 
% Apply f to another CHEBFUN3. 1-index contraction of two CHEBFUN3 objects
% will be a 4D function. 2-index contraction of two CHEBFUN3 objects will
% also be a CHEBFUN2. The first one is not possible in CHEBFUN and the
% second one is also not of interest. Here we check these. It should
% currently make an error:
g3D = chebfun3(@(x,y,z) x+y+z);
try
    h = f*g3D;
    pass(4) = false;
catch ME
    pass(4) = true;
end
 
end
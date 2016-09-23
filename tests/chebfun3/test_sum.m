function pass = test_sum(pref)
% Test chebfun3/sum.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e4*pref.cheb3Prefs.chebfun3eps;

ff = 'x.^2 + 4*y - 5*z';
dom = [11 14 7 10 5 9];
f = chebfun3(ff, dom);

exact1 = chebfun2(@(y,z) 12*y - 15*z + 471, dom(3:6));

pass(1)  = norm(sum(f) - exact1) < tol;

pass(2)  = norm(sum(f, [1]) - exact1) < tol;

exact2 = chebfun2(@(x,z) 3*x.^2 - 15*z + 102, [dom(1:2), dom(5:6)]);
pass(3)  = norm(sum(f, [2]) - exact2) < tol;

exact3 = chebfun2(@(x,y) 4*x.^2 + 16*y - 140, dom(1:4));
pass(4)  = norm(sum(f, [3]) - exact3) < tol;

% Make sure there is no bug when f has no x coordinate and we integrate
% over x:
dom = [-2 3, -2 2, -5 -3];
f = chebfun3(@(x,y,z) cos(y) + sin(z), dom);
g = sum(f);
exact = chebfun2(@(y,z) diff(dom(1:2)).*(cos(y) + sin(z)), dom(3:6));
pass(5)  = norm(g - exact) < tol;

% Make sure there is no bug when f has no y coordinate and we integrate
% over x:
f = chebfun3(@(x,y,z) cos(x) + sin(z), dom);
g = sum(f);
exact = chebfun2(@(y,z) sin(3)+sin(2) + 5.*sin(z), dom(3:6));
pass(6)  = norm(g - exact) < tol;

% Make sure there is no bug when f has no z coordinate and we integrate
% over x:
f = chebfun3(@(x,y,z) cos(x) + sin(y), dom);
g = sum(f);
exact = chebfun2(@(y,z) sin(3)+sin(2) + 5.*sin(y), dom(3:6));
pass(7)  = norm(g - exact) < tol;

% Make sure there is no bug when f has no x coordinate and we integrate over 
% y:
f = chebfun3(@(x,y,z) cos(y)+sin(z), dom);
g = sum(f, 2);
exact = chebfun2(@(x,z) 2*sin(2)+4.*sin(z), [dom(1:2), dom(5:6)]);
pass(8)  = norm(g - exact) < tol;

% Make sure there is no bug when f has no y coordinate and we integrate over 
% y:
f = chebfun3(@(x,y,z) cos(x)+sin(z), dom);
g = sum(f, 2);
exact = chebfun2(@(x,z) diff(dom(3:4)).*(cos(x)+sin(z)), [dom(1:2), dom(5:6)]);
pass(9)  = norm(g - exact) < tol;

% Make sure there is no bug when f has no z coordinate and we integrate over 
% y:
f = chebfun3(@(x,y,z) cos(x)+sin(y), dom);
g = sum(f, 2);
exact = chebfun2(@(x,z) 4*cos(x), [dom(1:2), dom(5:6)]);
pass(10)  = norm(g - exact) < tol;

% Make sure there is no bug when f has no x coordinate and we integrate over 
% z:
f = chebfun3(@(x,y,z) cos(y)+sin(z), dom);
g = sum(f, 3);
exact = chebfun2(@(x,y) 2*cos(y) + cos(5) - cos(3), dom(1:4));
pass(11)  = norm(g - exact) < tol;

% Make sure there is no bug when f has no y coordinate and we integrate over 
% z:
f = chebfun3(@(x,y,z) cos(x)+sin(z), dom);
g = sum(f, 3);
exact = chebfun2(@(x,y) 2*cos(x) + cos(5) - cos(3), dom(1:4));
pass(12)  = norm(g - exact) < tol;

% Make sure there is no bug when f has no z coordinate and we integrate over 
% z:
f = chebfun3(@(x,y,z) cos(x)+sin(y), dom);
g = sum(f, 3);
exact = chebfun2(@(x,y) diff(dom(5:6))*(cos(x)+sin(y)), dom(1:4));
pass(13)  = norm(g - exact) < tol;

end
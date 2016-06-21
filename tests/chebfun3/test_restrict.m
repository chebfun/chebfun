function pass = test_restrict( pref )
% This script checks the restriction of a Chebfun3 to smaller domains.

if ( nargin < 1 )
    pref = chebfunpref;
end 
tol = 1e3 * pref.cheb3Prefs.chebfun3eps;

ff = @(x,y,z) exp(x/2+y) + cos(x+z.^2);
dom = [-3 -2 -1 1 2 3];
f = chebfun3(ff, dom);

% Restricting a chebfun3 to a point.
val = f{-2.5,-2.5, 0,0, 2.5,2.5};  % should be ff(-2.5, 0, 2.5);
pass(1) = abs(val - ff(-2.5, 0, 2.5)) < tol;

% Restricting to a vertical line: 
f1D = f{dom(1),dom(2), 0,0, 2.5,2.5};
exact = chebfun(@(x) ff(x, 0,2.5), dom(1:2));
pass(2) = norm(f1D - exact) < tol;

% Restricting to a horizontal line: 
f1D = f{dom(1),dom(1), dom(3),dom(4), dom(5),dom(5)};  
exact = chebfun(@(y) ff(dom(1), y, dom(5)), dom(3:4));
pass(3) = norm(f1D - exact) < tol;

% Restricting to an oblique line:
f1D = f{dom(1),dom(1), dom(3),dom(3), dom(5),dom(6)};  
exact = chebfun(@(z) ff(dom(1), dom(3), z), dom(5:6));
pass(4) = norm(f1D - exact) < tol;

% Restricting a chebfun3 to a horizontal plane:
f2D = f{dom(1),dom(1), dom(3),dom(4), dom(5),dom(6)};  
exact = chebfun2(@(y,z) ff(dom(1), y, z), dom(3:6));
pass(5) = norm(f2D - exact) < tol;

% Restricting a chebfun3 to a lateral plane:
f2D = f{dom(1),dom(2), dom(3),dom(3), dom(5),dom(6)};  
exact = chebfun2(@(x,z) ff(x, dom(3), z), [dom(1:2) dom(5:6)]);
pass(6) = norm(f2D - exact) < tol;

% Restricting a chebfun3 to a frontal plane:
f2D = f{dom(1),dom(2), dom(3),dom(4), dom(5),dom(5)};  
exact = chebfun2(@(x,y) ff(x, y, dom(5)), dom(1:4));
pass(7) = norm(f2D - exact) < tol;

% Restricting a chebfun3 to a cuboid:
newDom = [-2.75 -2.25 -0.5 0.5 2.25 2.75];
g = f{newDom(1),newDom(2),newDom(3),newDom(4),newDom(5),newDom(6)};
exact = chebfun3(@(x,y,z) ff(x, y, z), newDom);
pass(8) = norm(g - exact) < tol;

% Restricting a chebfun3 to a cuboid that is not inside the original
% cuboid. This should give an error.
newDom = dom/2;
try
    g = f{newDom(1),newDom(2),newDom(3),newDom(4),newDom(5),newDom(6)};
    pass(9) = 0;
catch ME
    pass(9) = 1;
end

% Using subrefs for restricting
ref.type='{}'; 
ref.subs={-2.5,-2.5, 0,0, 2.5,2.5}; 
val = subsref(f, ref); 
exact = ff(-2.5, 0, 2.5);
pass(10) = abs(val - exact) < tol;

end
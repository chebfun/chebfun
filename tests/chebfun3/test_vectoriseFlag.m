function pass = test_vectoriseFlag(pref)
% Test the vectorise flag in the constructor. 

if ( nargin < 1 ) 
    pref = chebfunpref(); 
end 
tol = 10*pref.cheb3Prefs.chebfun3eps;

dom = [-1, 1, -1, 1, -1, 1];

% All these calls to the constructor should be the same: 
f1 = chebfun3(@(x,y,z) x); 
f2 = chebfun3(@(x,y,z) x, 'vectorize'); 
f3 = chebfun3(@(x,y,z) x, dom, 'vectorize');
f4 = chebfun3(@(x,y,z) x, 'vectorize', dom);
pass(1) = norm(f1 - f2) < tol;
pass(2) = norm(f1 - f3) < tol; 
pass(3) = norm(f1 - f4) < tol; 

% These should be the same: 
f1 = chebfun3(@(x,y,z) x.*y.*z); 
f2 = chebfun3(@(x,y,z) x*y*z, 'vectorize'); 
pass(4) = norm(f1 - f2) < tol;

end
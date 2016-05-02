function pass = test_sum(pref)
% Test chebfun3/sum.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 

tol = 1e4*pref.cheb3Prefs.chebfun3eps;
j = 1;

ff = 'x.^2 + 4*y - 5*z';
dom = [11 14 7 10 5 9];
f = chebfun3(ff, dom);

exact1 = chebfun2(@(y,z) 12*y - 15*z + 471, dom(3:6));

pass(j)  = norm(sum(f) - exact1) < tol;
j=j+1;

pass(j)  = norm(sum(f, [1]) - exact1) < tol;
j=j+1;

exact2 = chebfun2(@(x,z) 3*x.^2 - 15*z + 102, [dom(1:2), dom(5:6)]);
pass(j)  = norm(sum(f, [2]) - exact2) < tol;
j=j+1;

exact3 = chebfun2(@(x,y) 4*x.^2 + 16*y - 140, dom(1:4));
pass(j)  = norm(sum(f, [3]) - exact3) < tol;
j=j+1;

end
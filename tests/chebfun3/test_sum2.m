function pass = test_sum2(pref)
% Test chebfun3/sum2.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 

tol = 1e4*pref.cheb3Prefs.chebfun3eps;
j = 1;

ff = 'x.^2 + 4*y - 5*z';
dom = [11 14 7 10 5 9];
f = chebfun3(ff, dom);

exact12 = chebfun(@(z) 1719 - 45*z, dom(5:6));

pass(j)  = norm(sum2(f) - exact12) < tol;
j=j+1;

pass(j)  = norm(sum2(f, [1 2]) - exact12) < tol;
j=j+1;

pass(j)  = norm(sum2(f, [2 1]) - exact12) < tol;
j=j+1;

exact23 = chebfun(@(x) 12*x.^2 - 12, dom(1:2));
pass(j)  = norm(sum2(f, [2 3]) - exact23) < tol;
j=j+1;

pass(j)  = norm(sum2(f, [3 2]) - exact23) < tol;
j=j+1;

exact13 = chebfun(@(y) 48*y + 1464, dom(3:4));
pass(j)  = norm(sum2(f, [1 3]) - exact13) < tol;
j=j+1;

pass(j)  = norm(sum2(f, [3 1]) - exact13) < tol;
j=j+1;

end
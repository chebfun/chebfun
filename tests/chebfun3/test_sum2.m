function pass = test_sum2(pref)
% Test chebfun3/sum2.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 1e4*pref.cheb3Prefs.chebfun3eps;

ff = 'x.^2 + 4*y - 5*z';
dom = [11 14 7 10 5 9];
f = chebfun3(ff, dom);

exact12 = chebfun(@(z) 1719 - 45*z, dom(5:6));

pass(1)  = norm(sum2(f) - exact12) < tol;

pass(2)  = norm(sum2(f, [1 2]) - exact12) < tol;

pass(3)  = norm(sum2(f, [2 1]) - exact12) < tol;

exact23 = chebfun(@(x) 12*x.^2 - 12, dom(1:2));
pass(4)  = norm(sum2(f, [2 3]) - exact23) < tol;

pass(5)  = norm(sum2(f, [3 2]) - exact23) < tol;

exact13 = chebfun(@(y) 48*y + 1464, dom(3:4));
pass(6)  = norm(sum2(f, [1 3]) - exact13) < tol;

pass(7)  = norm(sum2(f, [3 1]) - exact13) < tol;

% A complex-valued function
f = chebfun3(@(x,y,z) 1i*x.^2 + 2i*y.^2 + z.^2);
exactZ = chebfun(@(z) 4*z.^2 + 4i);
pass(8)  = norm(sum2(f) - exactZ) < tol;

exactY = chebfun(@(y) 8i*y.^2 + 4/3*(1+1i));
pass(9)  = norm(sum2(f, [1 3]) - exactY) < tol;

exactX = chebfun(@(x) 1i*(4*x.^2+8/3) + 4/3);
pass(10)  = norm(sum2(f, [2 3]) - exactX) < tol;

end
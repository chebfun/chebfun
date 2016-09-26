function pass = test_zerofunction(pref)
% This test checks that a zero chebfun3 is being treated correctly for
% common commands. 

if ( nargin < 1 ) 
    pref = chebfunpref;
end
tol = 100 * pref.cheb3Prefs.chebfun3eps;

% construction
f = chebfun3(0); 
f = chebfun3(@(x,y,z) 0*x); 
f = chebfun3(0, [-2 2 -3 3 -4 4]);

% adding together a 0 chebfun3.
f = chebfun3(0); 
g = chebfun3(@(x,y,z) cos(x.*y.*z));
pass(1) = abs(norm(g + f) - norm(g)) < tol;

v = abs(f(pi/6, pi/6, pi/6)); 
v = v + rank(f);
pass(2) = (v==0);

v = abs(norm(sum(f)));
v = v + abs(sum3(f));

v = v + abs(norm(diff(f)));

v = v + abs(norm(diff(f,1,1)));

pass(3) = v==0;


% evaluation on an array. 
r = rand(10,8); 
v = f(r,r,r);
if all(size(v) == size(r))
    pass(4) = 1; 
end

end
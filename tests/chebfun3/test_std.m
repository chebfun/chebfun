function pass = test_std(pref)
% Test the chebfun3/max command. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 1e3*pref.cheb3Prefs.chebfun3eps;

% Note that std(chebfun(@(x) x.^2)) = sqrt(4/45).

f1 = chebfun3(@(x,y,z) x.^2 + y.*z); % f1 has x-std equal to sqrt(4/45).
f2 = chebfun3(@(x,y,z) y.^2 + x.*z); % f2 has y-std equal to sqrt(4/45). 
f3 = chebfun3(@(x,y,z) z.^2 + x.*y); % f3 has z-std equal to sqrt(4/45).

exact = chebfun2(@(x,y) sqrt(4/45)); % x-std of y.*z equals 0.

h1 = std(f1); 
h2 = std(f1,[],1); 
h3 = std(f2,[],2); 
h4 = std(f3,[],3); 

pass(1) = norm(h1 - exact) < tol;
pass(2) = norm(h2 - exact) < tol;
pass(3) = norm(h3 - exact) < tol;
pass(4) = norm(h4 - exact) < tol;

end
function pass = test_std2(pref)
% Test the chebfun3/max command. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 1e3*pref.cheb3Prefs.chebfun3eps;

% Note that std(chebfun(@(x) x.^2)) = sqrt(4/45).

f1 = chebfun3(@(x,y,z) x.^2 + z); % f1 has (x,y)-std equal to sqrt(4/45).
f2 = chebfun3(@(x,y,z) z.^2 + y); % f2 has (x,z)-std equal to sqrt(4/45). 
f3 = chebfun3(@(x,y,z) y.^2 + x); % f3 has (y,z)-std equal to sqrt(4/45).

exact = chebfun(@(x) sqrt(4/45)); % (x,y)-std of z equals 0.

h1 = std2(f1); 
h2 = std2(f1,[],[1,2]); 
h3 = std2(f2,[],[1,3]); 
h4 = std2(f3,[],[3, 2]); 

pass(1) = norm(h1 - exact) < tol;
pass(2) = norm(h2 - exact) < tol;
pass(3) = norm(h3 - exact) < tol;
pass(4) = norm(h4 - exact) < tol;

end
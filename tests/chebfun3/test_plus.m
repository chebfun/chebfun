function pass = test_plus(pref)
% Test chebfun3/plus.

if ( nargin < 1 )
    pref = chebfunpref; 
end 
tol = 1e5 * pref.cheb3Prefs.chebfun3eps;
j = 1;

dom = [-1 1 -1 1 -1 1; -2 2 -2 2 -2 2; -1 pi 0 2*pi -pi pi];
ff = @(x,y,z) cos(x.*y.*z);
gg = @(x,y,z) x + y + z + x.*y.*z;

for k = 1 : size(dom, 1)
    f = chebfun3(ff, dom(k, :));
    g = chebfun3(gg, dom(k, :));

    FplusG = chebfun3(@(x,y,z) ff(x,y,z) + gg(x, y, z), dom(k, :));
    
    tolk = norm(dom(k, :), inf) * tol;
    
    pass(j) = norm( f+g - FplusG ) < tolk;
    j = j + 1;
end

% Check if chebfun3/plus compresses the rank: 
f = chebfun3(@(x,y,z) x);
g = f + f; 
pass(4) = rank(g) == rank(f);

% Check if chebfun3/plus compresses the rank: 
f = chebfun3(@(x,y,z) cos(x.*y.*z)); 
g = f + f;
pass(5) = rank(g) <= rank(f);

% Check adding a function with a small vscale works: 
f = chebfun3(@(x,y,z) 1e-10*x);
g = f + f; 
pass(6) = rank(g) == rank(f);
pass(7) = abs(vscale(g) - 2e-10) < tol;

% Check adding a function with a large vscale works: 
f = chebfun3(@(x,y,z) 1e100*x); 
g = f + f;
pass(8) = rank(g) == rank(f);
pass(9) = abs(vscale(g) - 2e100)/2e100 < 2e-2;

% Does plus work fine with objects of different tech?
ff = @(x,y,z) sin(pi*x).*cos(pi*(x+y));
f = chebfun3(ff);
g = chebfun3(ff,'trig');
fmg = f - g;
pass(10) = norm(fmg) < tol;

end
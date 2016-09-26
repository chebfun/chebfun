function pass = test_syntax(pref)
% Check different syntax in CHEBFUN3V constructor.

if ( nargin < 1 )
    pref = chebfunpref;
end
tol = 1e5 * pref.cheb3Prefs.chebfun3eps;

dom = [-1 1 -1 1 -1 1; 
       -1 1 1 2 -5 -3];

for jj = 1 : size(dom, 1)    
    f = @(x,y,z) cos(x.*z) + sin(x.*y);
    g = @(x,y,z) cos(x.*y.*z);
    
    fcheb = chebfun3(f, dom(jj, :));
    gcheb = chebfun3(g, dom(jj, :));
    
    F1 = chebfun3v(f, g , dom(jj, :));
    F2 = chebfun3v({f; g}, dom(jj, :));
    F3 = chebfun3v(fcheb, gcheb);
    F4 = chebfun3v(fcheb, gcheb, dom(jj, :));
    
    pass(1, jj) = norm(F1 - F2) < tol; 
    pass(2, jj) = norm(F2 - F3) < tol;
    pass(3, jj) = norm(F3 - F4) < tol; 
end

end
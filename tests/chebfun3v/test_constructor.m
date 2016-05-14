function pass = test_constructor( pref )
% Test the Chebfun3v constructor when performing simple arithmetic
% operations.

if ( nargin < 1 )
    pref = chebfunpref;
end

tol = 100*pref.cheb3Prefs.chebfun3eps;

D = [-1 1 -1 1 -1 1; 
    -2 3 0 1 2 4];

for kk = 1:size(D, 1)
    dom = D(kk,:);
    % Check the constructor works with lots of different syntax:
    f1 = @(x,y,z) cos(x.*y.*z);
    f2 = @(x,y,z) sin(x+z) + cos(y);
    f3 = @(x,y,z) exp(x) - exp(-y) + z;
    g1 = chebfun3(f1, dom);
    g2 = chebfun3(f2, dom);
    g3 = chebfun3(f3, dom);
    
    % Two components:
    H1 = chebfun3v(f1, f2, dom);
    H2 = chebfun3v(f1, g2, dom);
    H3 = chebfun3v(g1, f2, dom);
    H4 = chebfun3v(g1, g2, dom);
    
    pass(1) = norm(H1 - H2) < tol;
    pass(2) = norm(H1 - H3) < tol;
    pass(3) = norm(H1 - H4) < tol;
    
    % Three components:
    H1 = chebfun3v(f1, f2, f3, dom);
    H2 = chebfun3v(f1, f2, g3, dom);
    H3 = chebfun3v(f1, g2, f3, dom);
    H4 = chebfun3v(g1, f2, f3, dom);
    H5 = chebfun3v(g1, g2, f3, dom);
    H6 = chebfun3v(g1, f2, g3, dom);
    H7 = chebfun3v(f1, g2, g3, dom);
    H8 = chebfun3v(g1, g2, g3, dom);
    
    pass(4) = norm(H1 - H2) < tol;
    pass(5) = norm(H1 - H3) < tol;
    pass(6) = norm(H1 - H4) < tol;
    pass(7) = norm(H1 - H5) < tol;
    pass(8) = norm(H1 - H6) < tol;
    pass(9) = norm(H1 - H7) < tol;
    pass(10) = norm(H1 - H8) < tol;
    
end

end
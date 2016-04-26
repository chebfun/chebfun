function pass = test_syntax( pref )
% Check the Chebfun2v constructor for different syntax.
% Alex Townsend, March 2013.

if ( nargin < 1 )
    pref = chebfunpref;
end
tol = 1e5 * pref.cheb2Prefs.chebfun2eps;

for jj = 1 : 2
    
    f = @(x,y,z) jj*sin(x.*y.*z);  % simple function.
    g = @(lam,th) exp((cos(lam).*sin(th)).^jj);
    h = @(x,y,z) cos(jj*x) + f(x,y,z);
    
    fsphere = spherefun(f);
    gsphere = spherefun(g);
    hsphere = spherefun(h);
    
    F1 = spherefunv( f, g, h );
    F2 = spherefunv( {f; g; h} );
    F3 = spherefunv( fsphere, gsphere, hsphere );
    
    pass(1, jj) = ( norm(F1 - F2) < tol ); 
    pass(2, jj) = ( norm(F2 - F3) < tol );
    
end
pass = pass(:)';

end
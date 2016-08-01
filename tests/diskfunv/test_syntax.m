function pass = test_syntax( pref )
% Check the Chebfun2v constructor for different syntax.
% Alex Townsend, March 2013.

if ( nargin < 1 )
    pref = chebfunpref;
end
tol = 1e5 * pref.cheb2Prefs.chebfun2eps;

for jj = 1 : 2
    
    f = @(x,y) jj*sin(x.*y);  % simple function.
    g = @(x,y ) jj*cos(x.*y);
    
    fd = diskfun(f);
    gd = diskfun(g);

    F1 = diskfunv( f, g);
    F2 = diskfunv( fd, gd);
    F3 =  [ fd; gd]; 
    
    pass(1, jj) = ( norm(F1 - F2) < tol ); 
    pass(2, jj) = ( norm(F2 - F3) < tol );
    
end
pass = pass(:)';

end
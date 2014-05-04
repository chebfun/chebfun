function pass = chebfun2v_syntax( pref )
% Check the Chebfun2v constructor for different syntax.
% Alex Townsend, March 2013.

if ( nargin < 1 )
    pref = chebfunpref;
end
tol = 1e5 * pref.cheb2Prefs.eps; 

j = 1;
D = [-1 1 -1 1; -1 1 1 2];


for jj = 1 : size(D, 1)
    
    f = @(x,y) cos(x) + sin(x.*y);  % simple function.
    g = @(x,y) cos(x.*y);
    
    fcheb = chebfun2(f, D(jj, :) );
    gcheb = chebfun2(g, D(jj, :) );
    
    F1 = chebfun2v( f, g , D(jj, :));
    F2 = chebfun2v( {f, g}, D(jj, :) );
    F3 = chebfun2v( fcheb, gcheb );
    F4 = chebfun2v( fcheb, gcheb, D(jj, :) );
    
    pass(j) = ( norm(F1 - F2) < tol ); j = j + 1; 
    pass(j) = ( norm(F2 - F3) < tol ); j = j + 1;
    pass(j) = ( norm(F3 - F4) < tol ); j = j + 1; 
    
end
end
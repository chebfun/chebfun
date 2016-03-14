function pass = test_plus( pref )
% Quick test for chebop2/plus. 

if ( nargin < 1 ) 
    pref = chebfunpref(); 
end
tol = pref.cheb2Prefs.chebfun2eps;

% Make a laplacian: 
N1 = chebop2(@(u) diff(u,2,2)); 
N2 = chebop2(@(x,y,u) diff(u,2,1));
N = N1 + N2; 

pass(1) = ( N.xorder == 2 );
pass(2) = ( N.yorder == 2 ); 
pass(3) = ( norm( cell2mat(N.coeffs) - [0 0 1; 0 0 0; 1 0 0] ) < tol);

end
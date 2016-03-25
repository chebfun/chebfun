function pass = test_laplacian( pref ) 
% Test LAPLACIAN 

if ( nargin == 0 ) 
    pref = chebfunpref; 
end

tol = 50*pref.cheb2Prefs.chebfun2eps;

% Check definition: 
F = chebfun2v(@(x,y) cos(x), @(x,y) sin(y)); 
lapF1 = diff(F(1), 2, 1) + diff(F(1), 2, 2);   
lapF2 = diff(F(2), 2, 1) + diff(F(2), 2, 2);  
lapF = laplacian(F); 
pass(1) = ( norm(lapF1 - lapF(1) ) < tol );
pass(2) = ( norm(lapF2 - lapF(2) ) < tol );


% Check definition: 
F = chebfun2v(@(x,y) cos(x), @(x,y) sin(y), @(x,y) x.*y + y.^2); 
lapF1 = diff(F(1), 2, 1) + diff(F(1), 2, 2);   
lapF2 = diff(F(2), 2, 1) + diff(F(2), 2, 2);
lapF3 = diff(F(3), 2, 1) + diff(F(3), 2, 2);
lapF = laplacian(F); 
pass(1) = ( norm(lapF1 - lapF(1) ) < tol );
pass(2) = ( norm(lapF2 - lapF(2) ) < tol );
pass(3) = ( norm(lapF3 - lapF(3) ) < tol );
end
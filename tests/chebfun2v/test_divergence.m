function pass = test_divergence( pref ) 
% Test DIVERGENCE

if ( nargin == 0 ) 
    pref = chebfunpref; 
end

tol = 100*pref.cheb2Prefs.chebfun2eps;

% Check definition: 
F = chebfun2v(@(x,y) cos(x), @(x,y) sin(y)); 
divF = diff(F(1), 1, 2) + diff(F(2), 1, 1);   
pass(1) = ( norm(divF - divergence(F) ) < tol );

% check third component is ignored: 
F = chebfun2v(@(x,y) cos(x), @(x,y) sin(y), @(x,y) cos(x)); 
divF = diff(F(1), 1, 2) + diff(F(2), 1, 1);   
pass(2) = ( norm(divF - divergence(F) ) < tol );

% Check divergence of gradient is laplacian: 
f = chebfun2(@(x,y) cos(x.*y)); 
pass(3) = ( norm(laplacian( f ) - divergence( gradient( f ) ) ) < tol );

end
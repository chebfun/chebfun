function pass = test_cross( pref ) 
% Test CROSS

if ( nargin == 0 ) 
    pref = chebfunpref; 
end

tol = 10*pref.eps; 

% Check definition: 
F = chebfun2v(@(x,y) cos(x), @(x,y) sin(y)); 
G = chebfun2v(@(x,y) x, @(x,y) y); 
crossF = F(1) .* G(2) - F(2) .* G(1);
pass(1) = ( norm(cross(F,G) - crossF) < tol );

% Check definition: 
F = chebfun2v(@(x,y) cos(x), @(x,y) sin(y), @(x,y) x.*y); 
G = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y ); 
crossF = [ F(2) .* G(3) - F(3) .* G(2) ; ...
          F(3) .* G(1) - F(1) .* G(3) ; ...
          F(1) .* G(2) - F(2) .* G(1) ];
pass(2) = ( norm(cross(F,G) - crossF) < tol );

end
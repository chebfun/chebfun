function pass = test_cross(pref)
% Test CROSS

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 10*pref.cheb3Prefs.chebfun3eps;

% Check definition: 
F = chebfun3v(@(x,y,z) cos(x), @(x,y,z) sin(y), @(x,y,z) exp(z)); 
G = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
crossF = [ F(2) .* G(3) - F(3) .* G(2) ; ...
          F(3) .* G(1) - F(1) .* G(3) ; ...
          F(1) .* G(2) - F(2) .* G(1) ];
pass(1) = norm(cross(F,G) - crossF) < tol;

% Check definition: 
F = chebfun3v(@(x,y,z) cos(x), @(x,y,z) sin(y), @(x,y,z) x.*y); 
G = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) x + y ); 
crossF = [ F(2) .* G(3) - F(3) .* G(2) ; ...
          F(3) .* G(1) - F(1) .* G(3) ; ...
          F(1) .* G(2) - F(2) .* G(1) ];
pass(2) = norm(cross(F,G) - crossF) < tol;

end
function pass = test_repeatedArithmetic( prefs )
% In version 4 repeated arithmetric lost quite a bit of accurate. 
% With PLUS based on QR, we are in better shape. Test this. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end
tol = 1e4*prefs.eps; 

% Add 50 functions together: 
f = chebfun2(@(x,y) cos(x.*y)); 
g = chebfun2(0);
for jj = 1:50 
    g = g + f; 
end
pass(1) = ( norm( g - 50*f ) < 10 * tol ); 

% Multiply 10 functions together: 
f = chebfun2(@(x,y) cos(x.*y)); 
g = f;
for jj = 1:10 
    g = g .* f; 
end
pass(2) = ( norm( g - f.^11 ) < tol );

% Multiply 20 functions together: 
f = chebfun2(@(x,y) sin(x.*y)); 
g = f;
for jj = 1:20 
    g = g .* f; 
end
pass(3) = ( norm( g - f.^21 ) < tol );

% Near cancellation not causing problems: 
f = chebfun2(@(x,y) cos(x.*y)); 
g = f + 1e-15*chebfun2(@(x,y) sin(x.*y.^2)); 
pass(4) = ( norm( f - g ) <  tol );

end
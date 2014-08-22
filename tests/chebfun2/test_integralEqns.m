function pass = test_integralEqns( pref ) 
% Test fred and volt

if ( nargin == 0) 
    pref = chebfunpref; 
end

tol = 1000*pref.eps; 
j = 1; 

% Fred on [-1 1 -1 1]
F = chebfun2(@(x,y) cos(x) + sin(y) ); 
v = chebfun(@(x) cos(x)); 
g = fred(F, v);
exact = chebfun(@(x) 2*sin(1)*cos(x))';
pass(j) = norm( g - exact ) < tol; j = j + 1; 

g = fred(F', v);
exact = chebfun(@(x) 2*sin(1)*sin(x) + 1 +sin(1)*cos(1))';
pass(j) = norm( g - exact ) < tol; j = j + 1; 


% Fred on [-3 4 -2 0]
F = chebfun2(@(x,y) cos(x) + sin(y), [-3 4 -2 0]); 
v = chebfun(@(x) cos(x), [-2 0]); 
g = fred(F, v);
exact = chebfun(@(x) sin(2)*cos(x)-sin(2).^2/2, [-3 4])';
pass(j) = norm( g - exact ) < tol; j = j + 1; 

v = chebfun(@(x) cos(x), [-3 4]); 
g = fred(F', v);
exact = chebfun(@(x) (sin(3) + sin(4))*sin(x) + .25*(14+sin(6)+sin(8)), [-2 0])';
pass(j) = norm( g - exact ) < tol; j = j + 1; 


% Volt on [-1 1 -1 1]
F = chebfun2(@(x,y) cos(x) + sin(y) ); 
v = chebfun(@(x) cos(x)); 
g = volt(F, v);
exact = chebfun(@(x) .25*(cos(2)-cos(2*x)) + (sin(x) + sin(1)).*cos(x))';
pass(j) = norm( g - exact ) < tol; j = j + 1; 

g = volt(F', v);
exact = chebfun(@(x) sin(x).*(sin(x) + sin(1)) + .25*(2*x+sin(2*x) + 2 + sin(2)))';
pass(j) = norm( g - exact ) < tol; j = j + 1; 

% Test fred construction with the deprecated style.
K = @(s,t) exp(-abs(s-t));
F = fred( K, domain([-1,1]) );
pass(j) = all(F.domain == [-1, 1]); j = j + 1;
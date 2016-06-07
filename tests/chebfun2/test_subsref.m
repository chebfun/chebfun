function pass = test_subsref( pref )

if ( nargin < 1 ) 
    pref = chebfunpref; 
end

tol = 1000*pref.cheb2Prefs.chebfun2eps;

% Evaluation 
f = chebfun2(@(x,y) sin(x.*(y-.1)), [-2 2 -3 4]);

pass(1) = ( norm( f( pi/4, pi/6 ) - sin(pi/4.*(pi/6-.1)) ) < tol ); 

slice = chebfun(@(x) sin(x.*(pi/6-.1)), [-2 2] );
slice2 = chebfun(@(x) sin(x.*(pi/4-.1)), [-2 2] ); 
pass(2) = ( norm( f( :, [pi/6 pi/4] ) - [slice slice2]' ) < tol );

slice = chebfun(@(y) sin(pi/4.*(y-.1)), [-3 4] ); 
slice2 = chebfun(@(y) sin(pi/6.*(y-.1)), [-3 4] ); 
pass(3) = ( norm( f( pi/4, : ) - slice ) < tol ); 
pass(4) = ( norm( f( [pi/4 pi/6], : ) - [slice slice2] ) < tol ); 

pass(5) = ( norm( f( :, : ) - f ) < tol );

% Test evaluation syntax for chebfun inputs. 
f = chebfun2(@(x,y) x.*y); 
c1 = chebfun(@(t) 1 + 0*t);
c2 = chebfun(@(t) -.3 + 0*t);
pass(6) = ( norm( f(c1, c2) +.3 ) < tol );
pass(7) = ( norm( f(c1 + 1i*c2) +.3 ) < tol );
pass(8) = ( norm( feval(f, c1, c2) - f(c1, c2) ) < tol );

% GET properties 
f = chebfun2(@(x,y) x);  
pass(9) = ( norm( f.rows - chebfun(@(x) x) ) < tol ); 
pass(10) = ( norm( f.domain - [-1 1 -1 1] ) < tol ); 

% Restriction 
f = chebfun2(@(x,y) sin(x.*(y-.1)), [-2 2 -3 4]);
g = f{-1,1,-.5,.25}; 
exact = chebfun2(@(x,y) sin(x.*(y-.1)), [-1,1,-.5,.25]);
pass(11) = ( norm( g -  exact ) < 10*tol );

end
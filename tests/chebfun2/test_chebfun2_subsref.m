function pass = test_chebfun2_subsref( pref )

if ( nargin < 1 ) 
    pref = chebfunpref; 
end

tol = 1000*pref.cheb2Prefs.eps; 
j = 1; 

% Evaluation 
f = chebfun2(@(x,y) sin(x.*(y-.1)), [-2 2 -3 4]);

pass(j) = ( norm( f( pi/4, pi/6 ) - sin(pi/4.*(pi/6-.1)) ) < tol ); j = j + 1; 

slice = chebfun(@(x) sin(x.*(pi/6-.1)), [-2 2] );
slice2 = chebfun(@(x) sin(x.*(pi/4-.1)), [-2 2] ); 
pass(j) = ( norm( f( :, [pi/6 pi/4] ) - [slice slice2]' ) < tol ); j = j + 1;

slice = chebfun(@(y) sin(pi/4.*(y-.1)), [-3 4] ); 
slice2 = chebfun(@(y) sin(pi/6.*(y-.1)), [-3 4] ); 
pass(j) = ( norm( f( pi/4, : ) - slice ) < tol ); j = j + 1;
pass(j) = ( norm( f( [pi/4 pi/6], : ) - [slice slice2] ) < tol ); j = j + 1;

pass(j) = ( norm( f( :, : ) - f ) < tol ); j = j + 1;

% Test evaluation syntax for chebfun inputs. 
f = chebfun2(@(x,y) x.*y); 
c1 = chebfun(@(t) 1 + 0*t);
c2 = chebfun(@(t) -.3 + 0*t);
pass(j) = ( norm( f(c1, c2) +.3 ) < tol );
pass(j) = ( norm( f(c1 + 1i*c2) +.3 ) < tol );
pass(j) = ( norm( feval(f, c1, c2) - f(c1, c2) ) < tol );

% GET properties 
f = chebfun2(@(x,y) x); 
pass(j) = ( norm( abs(f.cols) - 1 ) < tol ); j = j + 1; 
pass(j) = ( norm( f.rows - chebfun(@(x) x) ) < tol ); j = j + 1; 
pass(j) = ( norm( f.domain - [-1 1 -1 1] ) < tol ); j = j + 1; 
pass(j) = ( norm( abs(f.pivotValues) - 1 ) < tol ); j = j + 1; 

% Restriction 
f = chebfun2(@(x,y) sin(x.*(y-.1)), [-2 2 -3 4]);
g = f{-1,1,-.5,.25}; 
exact = chebfun2(@(x,y) sin(x.*(y-.1)), [-1,1,-.5,.25]);
pass(j) = ( norm( g -  exact ) < 10*tol ); j = j + 1;

end
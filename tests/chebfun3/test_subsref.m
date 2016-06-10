function pass = test_subsref(pref)

if ( nargin < 1 ) 
    pref = chebfunpref; 
end

tol = 1e4*pref.cheb3Prefs.chebfun3eps;
j = 1;

% Evaluation 
ff = @(x,y,z) sin(x.*(y-.1).*(z+.2));
dom = [-2 2 -3 4 -pi pi];
f = chebfun3(ff, dom);
pass(j) = abs(f(pi/4, pi/6, pi/3) - ff(pi/4, pi/6, pi/3)) < tol;
j = j+1;

cross1 = chebfun2(@(x,y) sin(x.*(y-.1).*(pi/6+.2)), dom(1:4));
cross2 = chebfun2(@(x,y) sin(x.*(y-.1).*(pi/4+.2)), dom(1:4));
pass(j) = norm(f(:, :, pi/6) - cross1) < tol && ...
            norm(f( :, :, pi/4) - cross2) < tol;
j = j+1;

cross1 = chebfun2(@(y,z) sin(pi/4.*(y-.1).*(z+.2)), dom(3:6)); 
cross2 = chebfun2(@(y,z) sin(pi/6.*(y-.1).*(z+.2)), dom(3:6)); 
pass(j) = norm(f(pi/4, :, :) - cross1) < tol; 
j = j+1;
pass(j) = norm(f(pi/6, :, :) - cross2) < tol; 
j= j+1;

pass(j) = norm( f(:, :, :) - f) < tol;
j = j+1;

% Test evaluation syntax for chebfun inputs. 
f = chebfun3(@(x,y,z) x.*y.*z); 
c1 = chebfun(@(t) 1 + 0*t);
c2 = chebfun(@(t) -.3 + 0*t);
c3 = chebfun(@(t) 0.5 + 0*t);
pass(j) = norm(f(c1, c2, c3) - 1*(-0.3)*0.5 ) < tol;
j = j+1;
pass(j) = norm(feval(f, c1, c2, c3) - f(c1, c2, c3)) < tol;
j = j+1;

% GET properties 
f = chebfun3(@(x,y,z) x);  
pass(j) = norm(f.cols - chebfun(@(x) x) ) < tol;
j = j+1;
pass(j) = norm(f.domain - [-1 1 -1 1 -1 1] ) < tol;
j = j+1;

end
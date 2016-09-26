function pass = test_subsref( pref )

if ( nargin < 1 ) 
    pref = chebfunpref; 
end

tol = 1000*pref.cheb2Prefs.chebfun2eps;

% Evaluation at a point and test flags
ff = @(x,y) sin(x.*(y-.1));
f = diskfun(ff);
x = 1/sqrt(3); y = 1/sqrt(3);
pass(1) = ( abs( f(x,y) - ff(x,y) ) < tol ); 
% try polar coords
r = 1/3; t = pi/7; 
pass(2) = (abs(f(t,r, 'polar')-ff(r.*cos(t), r.*sin(t))) < tol); 
pass(3) = ( abs( f(x,y, 'cart') - ff(x,y) ) < tol ); 

% Slices in r
ff = @(t,r) exp(r.^2.*sin(t).*cos(t));
f = diskfun(ff, 'polar');
th1 = 0.2;
slice1 = chebfun(@(rad) feval(ff,th1, rad));
th2 = 1.7;
slice2 = chebfun(@(rad) feval(ff,th2, rad));
pass(4) = ( norm( [f(th1,:) f(th2,:)] - [slice1 slice2] ) < tol );
pass(5) = ( norm( f([th1 th2],:) - [slice1 slice2] ) < tol );
% Slices in theta
r1 = 0.2;
slice1 = chebfun(@(th) feval(ff,th, r1), [-pi, pi], 'trig');
r2 = .7;
slice2 = chebfun(@(th) feval(ff,th, r2), [-pi, pi], 'trig');
pass(6) = ( norm( f(:,[r1 r2]) - [slice1 slice2] ) < tol );
%choose negative r
r3 = -.7; 
pass(7) = (norm(f(:, r2)-f(:,r3))<tol);


% GET properties 
f = diskfun(@(t,r) r.*cos(t), 'polar');  
pass(8) = ( norm(f.rows - chebfun(@(t) cos(t),[-pi pi],'trig')) < tol || ...
    norm(f.rows + chebfun(@(t) cos(t),[-pi pi],'trig')) < tol ); 
pass(9) = ( norm(f.cols - chebfun(@(r) r)) < tol || ...
    norm(f.cols + chebfun(@(r) r)) < tol ); 
pass(10) = ( norm( f.domain - [-pi pi 0 1] ) < tol ); 

end
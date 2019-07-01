function pass = test_subsref( pref )

if ( nargin < 1 ) 
    pref = chebfunpref; 
end

tol = 1000*pref.cheb2Prefs.chebfun2eps;

% Evaluation at x,y,z
ff = @(x,y,z) z.*sin(x.*(y-.1));
f = spherefun(ff);
x = 1/sqrt(3); y = 1/sqrt(3); z = -1/sqrt(3);
pass(1) = ( abs( f(x,y,z) - ff(x,y,z) ) < tol ); 

% Evaluation at (lambda,theta)
ff = @(lam,th) exp(cos(lam).*sin(th).*cos(th));
f = spherefun(ff);
lam = pi/3; th = pi/4;
pass(2) = ( abs( f(lam,th ) - ff(lam,th) ) < tol ); 

% Slices in theta
ff = @(lam,th) exp(cos(lam).*sin(th).*cos(th));
f = spherefun(ff);
th1 = 0.2;
slice1 = chebfun(@(lam) feval(ff,lam,th1),[-pi pi], 'trig');
th2 = 0.7;
slice2 = chebfun(@(lam) feval(ff,lam,th2),[-pi pi], 'trig');
pass(3) = ( norm( f(:,[th1 th2]).' - [slice1 slice2] ) < tol );

% Slices in lambda
ff = @(lam,th) exp(cos(lam).*sin(th).*cos(th));
f = spherefun(ff);
lam1 = 0.2;
slice1 = chebfun(@(th) feval(ff,lam1,th),[-pi pi], 'trig');
lam2 = 0.7;
slice2 = chebfun(@(th) feval(ff,lam2,th),[-pi pi], 'trig');
pass(4) = ( norm( f([th1 th2],:) - [slice1 slice2] ) < tol );

% Slice in z
ff = @(x,y,z) z.*sin(x.*(y-.1));
f = spherefun(ff);
z = 0.1; th = acos(z);
slice = chebfun(@(lam) feval(ff,cos(lam).*sin(th),sin(lam).*sin(th),z),[-pi pi], 'trig');
pass(5) = ( norm( f(:,:,z) - slice ) < tol ); 

% Slice in x
ff = @(x,y,z) z.*sin(x.*(y-.1));
f = spherefun(ff);
x = 0.1;
slice = chebfun(@(t) feval(ff,x,sqrt(1-x.^2).*cos(t),sqrt(1-x.^2).*sin(t)),[-pi pi], 'trig');
pass(6) = ( norm( f(x,:,:) - slice ) < tol ); 

% Slice in y
ff = @(x,y,z) z.*sin(x.*(y-.1));
f = spherefun(ff);
y = -0.4;
slice = chebfun(@(t) feval(ff,sqrt(1-y.^2).*cos(t),y,sqrt(1-y.^2).*sin(t)),[-pi pi], 'trig');
pass(7) = ( norm( f(:,y,:) - slice ) < tol ); 

% GET properties 
f = spherefun(@(lam,th) cos(lam).*sin(th));  
pass(8) = ( norm(f.rows - chebfun(@(lam) cos(lam),[-pi pi],'trig')) < tol || ...
    norm(f.rows + chebfun(@(lam) cos(lam),[-pi pi],'trig')) < tol ); 
pass(9) = ( norm(f.cols - chebfun(@(th) sin(th),[-pi pi],'trig')) < tol || ...
    norm(f.cols + chebfun(@(th) sin(th),[-pi pi],'trig')) < tol ); 
pass(10) = ( norm( f.domain - [-pi pi 0 pi] ) < tol ); 

% Composition of a spherefun with a chebfun (1 and 3 columns):
f = spherefun(@(x,y,z) z + sin(pi*x.*y));
g = chebfun(@(t) t.^2, [ -1.5, 1.5 ]);
h_true = spherefun(@(x,y,z) (z + sin(pi*x.*y)).^2);
h = g(f);
pass(11) = ( norm(h - h_true) < tol );

G = chebfun(@(t) [ t.^2, t, -t.^2 ], [ -1.5, 1.5 ]);
H_true = spherefunv(@(x,y,z) (z + sin(pi*x.*y)).^2, ...
    @(x,y,z) z + sin(pi*x.*y), @(x,y,z) -(z + sin(pi*x.*y)).^2);
H = G(f);
pass(12) = ( norm(H - H_true) < tol );

end
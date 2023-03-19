function pass = test_subsref(pref)

if ( nargin < 1 )
    pref = chebfunpref;
end
tol = 1e4*pref.cheb3Prefs.chebfun3eps;

% Evaluation
ff = @(x,y,z) sin(x.*(y-.1).*(z+.2));
dom = [-2 2 -3 4 -pi pi];
f = chebfun3(ff, dom);
pass(1) = abs(f(pi/4, pi/6, pi/3) - ff(pi/4, pi/6, pi/3)) < tol;

cross1 = chebfun2(@(x,y) sin(x.*(y-.1).*(pi/6+.2)), dom(1:4));
cross2 = chebfun2(@(x,y) sin(x.*(y-.1).*(pi/4+.2)), dom(1:4));
pass(2) = norm(f(:, :, pi/6) - cross1) < tol && ...
    norm(f( :, :, pi/4) - cross2) < tol;

cross1 = chebfun2(@(y,z) sin(pi/4.*(y-.1).*(z+.2)), dom(3:6));
cross2 = chebfun2(@(y,z) sin(pi/6.*(y-.1).*(z+.2)), dom(3:6));
pass(3) = norm(f(pi/4, :, :) - cross1) < tol;

pass(4) = norm(f(pi/6, :, :) - cross2) < tol;

pass(5) = norm(f(:, :, :) - f) < tol;

% Test evaluation syntax for chebfun inputs.
f = chebfun3(@(x,y,z) x.*y.*z);
c1 = chebfun(@(t) 1 + 0*t);
c2 = chebfun(@(t) -.3 + 0*t);
c3 = chebfun(@(t) 0.5 + 0*t);
pass(6) = norm(f(c1, c2, c3) - 1*(-0.3)*0.5 ) < tol;

pass(7) = norm(feval(f, c1, c2, c3) - f(c1, c2, c3)) < tol;

% GET properties
f = chebfun3(@(x,y,z) x);
pass(8) = norm((f.core .* f.rows .* f.tubes)*f.cols - chebfun(@(x) x) ) < tol; 

pass(9) = norm(f.domain - [-1 1 -1 1 -1 1] ) < tol;

% Test evaluations at wrong type of inputs. The following two should give
% an error:
f = chebfun3(@(x,y,z) x+y+z);
try
    v = f(0.5);
    pass(10) = false;
catch ME
    pass(10) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN3:subsref:inputs');
end

try
    v = f(0.5, 0.5);
    pass(11) = false;
catch ME
    pass(11) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN3:subsref:inputs');
end

% Test composition with a CHEBFUN2V or three CHEBFUN2
g = chebfun3(@(x,y,z) x + y + z, [ -1, 1, -1, 1, -2, 2 ]);
F = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y);
h = g(F);
h_true = chebfun2(@(x,y) 2*x + 2*y);
pass(12) = ( norm(h - h_true) < tol );

g = chebfun3(@(x,y,z) x + y + z, [ -1, 1, -1, 1, -2, 2 ]);
f1 = chebfun2(@(x,y) x);
f2 = chebfun2(@(x,y) y);
f3 = f1 + f2;
h = g(f1, f2, f3);
h_true = chebfun2(@(x,y) 2*x + 2*y);
pass(13) = ( norm(h - h_true) < tol );

% Test composition with a CHEBFUN3V or three CHEBFUN3:
g = chebfun3(@(x,y,z) x + y + z);
F = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
h = g(F);
pass(14) = ( norm(h - g) < tol );

g = chebfun3(@(x,y,z) x+y+z);
f1 = chebfun3(@(x,y,z) x);
f2 = chebfun3(@(x,y,z) y);
f3 = chebfun3(@(x,y,z) z);
h = g(f1, f2, f3);
pass(15) = ( norm(h - g) < tol );

% Test composition with one inf by 3 CHEBFUN:
F = chebfun(@(t) [t, t, t]);
g = chebfun3(@(x,y,z) x + y + z);
h = g(F);
h_true = chebfun(@(t) 3*t);
pass(16) = ( norm(h - h_true) < tol );

% Test composition with three CHEBFUNs:
f = chebfun(@(t) t);
g = chebfun3(@(x,y,z) x + y + z);
h = g(f, f, f);
h_true = chebfun(@(t) 3*t);
pass(17) = ( norm(h - h_true) < tol );

% Test composition with a SPHEREFUNV:
f = spherefunv(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
g = chebfun3(@(x,y,z) x.^2 + y.^2 + z.^2);
h = g(f);
h_true = spherefun(@(x,y,z) 0*x + 1);
pass(18) = ( norm(h - h_true) < tol );

end
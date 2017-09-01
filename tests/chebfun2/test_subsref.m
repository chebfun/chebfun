function pass = test_subsref(pref)
% Test SUBSREF()

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
pass(8) = ( norm( compose([ c1, c2 ], f) - f(c1, c2) ) < tol );

% GET properties
f = chebfun2(@(x,y) x);
pass(9) = ( norm( f.rows - chebfun(@(x) x) ) < tol );
pass(10) = ( norm( f.domain - [-1 1 -1 1] ) < tol );

% Restriction
f = chebfun2(@(x,y) sin(x.*(y-.1)), [-2 2 -3 4]);
g = f{-1,1,-.5,.25};
exact = chebfun2(@(x,y) sin(x.*(y-.1)), [-1,1,-.5,.25]);
pass(11) = ( norm( g -  exact ) < 10*tol );

% Evaluation of complex-valued Chebfun2 objects, #1956
f = chebfun2(@(z) z);
pass(12) = norm(f(1i) -  1i) < tol; 
pass(13) = norm(f(1) -  1) < tol; 

% Test chebfun2( chebfun2v ) composition.
f = chebfun2(@(x,y) cos(x.*y));
F = chebfun2v(@(x,y) x, @(x,y) y);
pass(14) = ( norm(f(F) - f) < tol );

f = chebfun2(@(x,y) cos(x.*y), [-3,4,-2,6]);
g = chebfun2(@(x,y) cos(x.*y), [-2,2,-2,2]);
F = chebfun2v(@(x,y) x, @(x,y) y, [-2,2,-2,2]);
pass(15) = ( norm(f(F) - g) < tol );

% Composition of a CHEBFUN2 with a CHEBFUN3V.
F = chebfun3v(@(x,y,z) x, @(x,y,z) y);
g = chebfun2(@(x,y) x + y);
h = g(F);
h_true = chebfun3(@(x,y,z) x + y);
pass(16) = ( norm(h - h_true) < tol );

% Composition g(f) of a CHEBFUN2 g with a complex CHEBFUN2 f, interpeted as
% [real(f); imag(f)].
f = chebfun2(@(x,y) x + 1i*y);
g = chebfun2(@(x,y) x + y);
h = g(f);
h_true = g;
pass(17) = ( norm(h - h_true) < tol );

% Composition g(f) of a CHEBFUN2 g with a complex CHEBFUN3 f, interpeted as
% [real(f); imag(f)].
f = chebfun3(@(x,y,z) x + 1i*y);
g = chebfun2(@(x,y) x + y);
h = g(f);
h_true = chebfun3(@(x,y,z) x + y);
pass(18) = ( norm(h - h_true) < tol );

% Composition g(f) of a CHEBFUN2 g with a real CHEBFUN3 f, interpeted as
% [real(f); imag(f)].
f = chebfun3(@(x,y,z) x);
g = chebfun2(@(x,y) x + y);
h = g(f);
h_true = chebfun3(@(x,y,z) x);
pass(19) = ( norm(h - h_true) < tol );

% Composition g(f1, f2) of a CHEBFUN2 g and two CHEBFUN3.
f1 = chebfun3(@(x,y,z) x);
f2 = chebfun3(@(x,y,z) y);
g = chebfun2(@(x,y) x + y);
h = g(f1, f2);
h_true = chebfun3(@(x,y,z) x + y);
pass(20) = ( norm(h - h_true) < tol );

% Composition g(f1, f2) of a CHEBFUN2 g and two CHEBFUN2.
f1 = chebfun2(@(x,y) x);
f2 = chebfun2(@(x,y) y);
g = chebfun2(@(x,y) x + y);
h = g(f1, f2);
h_true = chebfun2(@(x,y) x + y);
pass(21) = ( norm(h - h_true) < tol );

% Test composition with one inf by 2 CHEBFUN:
F = chebfun(@(t) [t, t]);
g = chebfun2(@(x,y) x + y);
h = g(F);
h_true = chebfun(@(t) 2*t);
pass(22) = ( norm(h - h_true) < tol );

% Test composition with two CHEBFUNs:
f = chebfun(@(t) t);
g = chebfun2(@(x,y) x + y);
h = g(f, f);
h_true = chebfun(@(t) 2*t);
pass(23) = ( norm(h - h_true) < tol );

% Test composition with a complex-valued CHEBFUN:
f = chebfun(@(t) t + 2i*t);
g = chebfun2(@(x,y) x + y);
h = g(f);
h_true = chebfun(@(t) 3*t);
pass(24) = ( norm(h - h_true) < tol );

% Test composition with a SPHEREFUN:
f = spherefun(@(x,y,z) x + y + z);
g = chebfun2(@(x,y) x + y, [-2, 2, -2, 2]);
h = g(f);
h_true = f;
pass(25) = ( norm(h - h_true) < tol );

% Test composition with a DISKFUN:
f = diskfun(@(x,y) x + y);
g = chebfun2(@(x,y) x.^2 + y, [-2, 2, -2, 2]);
h = g(f);
h_true = diskfun(@(x,y) (x + y).^2);
pass(26) = ( norm(h - h_true) < tol );

% Test composition with a DISKFUNV:
F = diskfunv(@(x,y) x, @(x,y) y);
g = chebfun2(@(x,y) x + y);
h = g(F);
h_true = diskfun(@(x,y) x + y);
pass(27) = ( norm(h - h_true) < tol );

end
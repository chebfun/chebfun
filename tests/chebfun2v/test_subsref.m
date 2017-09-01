function pass = test_subsref( pref )
% Test Chebfun2v subsref() command.

if ( nargin == 0)
    pref = chebfunpref;
end
tol = 1000*pref.cheb2Prefs.chebfun2eps;

% Test recursive subsref:
f = chebfun2(@(x, y) sin(x.*y));
F = [f ; f ; f];

G = F(1);
exact = G.pivotValues;

pass(1) = norm( exact - F(1).pivotValues ) < tol;

% Composition of a CHEBFUN2V (2 components) with two CHEBFUN2 objects.
f1 = chebfun2(@(x,y) x + 1);
f2 = chebfun2(@(x,y) 2*y);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y, [ 0, 2, -2, 2 ]);
H = G(f1, f2);
H_true = chebfun2v(@(x,y) x + 2*y + 1, @(x,y) x - 2*y + 1);
pass(2) = ( norm(H - H_true) < tol );
pass(3) = ( norm(H(1, 1) - [4; 0]) < tol );

% Composition of a CHEBFUN2V (3 components) with two CHEBFUN2 objects.
f1 = chebfun2(@(x,y) x + 1);
f2 = chebfun2(@(x,y) 2*y);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y, @(x,y) x.*y , [ 0, 2, -2, 2 ]);
H = G(f1, f2);
H_true = chebfun2v(@(x,y) x + 2*y + 1, @(x,y) x - 2*y + 1, ...
    @(x,y) 2*(x + 1).*y);
pass(4) = ( norm(H - H_true) < tol );

% Composition of a CHEBFUN2V with two CHEBFUN3 objects.
f1 = chebfun3(@(x,y,z) x);
f2 = chebfun3(@(x,y,z) y + z);
G = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y, [-1, 1, -2, 2 ]);
H = G(f1, f2);
H_true = chebfun3v(@(x,y,z) x, @(x,y,z) y + z, @(x,y,z) x + y + z);
pass(5) = ( norm(H - H_true) < tol );

% Composition of a CHEBFUN2V (2 components) with a CHEBFUN2V.
F = chebfun2v(@(x,y) x, @(x,y) y);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y);
H = G(F);
pass(6) = ( norm(H - G) < tol );

% Composition of a CHEBFUN2V (3 components) with a CHEBFUN2V.
F = chebfun2v(@(x,y) x, @(x,y) y);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y, @(x,y) x.*y);
H = G(F);
pass(7) = ( norm(H - G) < tol );

% Composition of a CHEBFUN2V (2 components) with a CHEBFUN3V.
F = chebfun3v(@(x,y,z) x, @(x,y,z) y + z);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y, [ -1, 1, -2, 2 ]);
H = G(F);
H_true = chebfun3v(@(x,y,z) x + y + z, @(x,y,z) x - y - z);
pass(8) = ( norm(H - H_true) < tol );

% Composition of a CHEBFUN2V (3 components) with a CHEBFUN3V.
F = chebfun3v(@(x,y,z) x, @(x,y,z) y + z);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y, @(x,y) x.*y, [ -1, 1, -2, 2 ]);
H = G(F);
H_true = chebfun3v(@(x,y,z) x + y + z, @(x,y,z) x - y - z, @(x,y,z) x.*(y + z));
pass(9) = ( norm(H - H_true) < tol );

% Composition of a CHEBFUN2V with a CHEBFUN2, which will be interpreted as the
% CHEBFUN2V [real(f); imag(f)].  (Regardless if f is real or complex.)
f = chebfun2(@(x,y) x + 1i*y);
G = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y);
H = G(f);
pass(10) = ( norm(H - G) < tol );

% With real f:
f = chebfun2(@(x,y) x);
G = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y);
H = G(f);
H_true = chebfun2v(@(x,y) x, @(x,y) 0*x, @(x,y) x);
pass(11) = ( norm(H - H_true) < tol );

% Composition of a CHEBFUN2V with a CHEBFUN3, which will be interpreted as the
% CHEBFUN3V [real(f); imag(f)].  (Regardless if f is real or complex.)
f = chebfun3(@(x,y,z) x + y + 1i*z);
G = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y, [-2, 2, -1, 1]);
H = G(f);
H_true = chebfun3v(@(x,y,z) x + y, @(x,y,z) z, @(x,y,z) x + y + z);
pass(12) = ( norm(H - H_true) < tol );

% With real f:
f = chebfun3(@(x,y,z) x + y);
G = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y, [ -2, 2, -2, 2 ]);
H = G(f);
H_true = chebfun3v(@(x,y,z) x + y, @(x,y,z) 0*y, @(x,y,z) x + y);
pass(13) = ( norm(H - H_true) < tol );

% Test composition with one inf by 2 CHEBFUN:
F = chebfun(@(t) [ t, t ]);
G = chebfun2v(@(x,y) x + y, @(x,y) x);
H = G(F);
t = chebfun(@(t) t);
H_true = [ 2*t, t ];
pass(14) = ( norm(H - H_true) < tol );

% Test composition with two CHEBFUNs:
f = chebfun(@(t) t);
G = chebfun2v(@(x,y) x + y, @(x,y) x);
H = G(f, f);
H_true = [ 2*f, f ];
pass(15) = ( norm(H - H_true) < tol );

% Test composition with a complex-valued CHEBFUN:
f = chebfun(@(t) t + 1i*t);
G = chebfun2v(@(x,y) x + y, @(x,y) x);
H = G(f);
t = chebfun(@(t) t);
H_true = [ 2*t, t ];
pass(16) = ( norm(H - H_true) < tol );

% Test composition with a SPHEREFUN:
f = spherefun(@(x,y,z) x);
G = chebfun2v(@(x,y) x + y, @(x,y) 2*x, @(x,y) y.^2);
H = G(f);
H_true = spherefunv(@(x,y,z) x, @(x,y,z) 2*x, @(x,y,z) 0*x);
pass(17) = ( norm(H - H_true) < tol );

% Test composition with a DISKFUN:
f = diskfun(@(x,y) x + y);
G = chebfun2v(@(x,y) x, @(x,y) y, [-2, 2, -2, 2]);
H = G(f);
H_true = diskfunv(@(x,y) x + y, @(x,y) 0*x);
pass(18) = ( norm(H - H_true) < tol );

% Test composition with a DISKFUNV:
F = diskfunv(@(x,y) x, @(x,y) y);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y);
H = G(F);
H_true = diskfunv(@(x,y) x + y, @(x,y) x - y);
pass(19) = ( norm(H - H_true) < tol );

end
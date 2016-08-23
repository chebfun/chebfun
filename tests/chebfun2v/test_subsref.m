function pass = test_subsref( pref )
% Test Chebfun2v subsref() command.

if ( nargin == 0)
    pref = chebfunpref;
end

tol = 1000*pref.cheb2Prefs.chebfun2eps;

j = 1;

% Test recursive subsref:
f = chebfun2(@(x, y) sin(x.*y));
F = [f ; f ; f];

G = F(1);
exact = G.pivotValues;

pass(j) = norm( exact - F(1).pivotValues ) < tol;
j = j+1;

% Composition of a CHEBFUN2V (2 components) with two CHEBFUN2 objects.
f1 = chebfun2(@(x,y) x + 1);
f2 = chebfun2(@(x,y) 2*y);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y, [ 0, 2, -2, 2 ]);
H = G(f1, f2);
H_true = chebfun2v(@(x,y) x + 2*y + 1, @(x,y) x - 2*y + 1);
pass(j) = ( norm(H - H_true) < tol );
j = j+1;
pass(j) = ( norm(H(1, 1) - [4; 0]) < tol );
j = j+1;

% Composition of a CHEBFUN2V (3 components) with two CHEBFUN2 objects.
f1 = chebfun2(@(x,y) x + 1);
f2 = chebfun2(@(x,y) 2*y);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y, @(x,y) x.*y , [ 0, 2, -2, 2 ]);
H = G(f1, f2);
H_true = chebfun2v(@(x,y) x + 2*y + 1, @(x,y) x - 2*y + 1, ...
    @(x,y) 2*(x + 1).*y);
pass(j) = ( norm(H - H_true) < tol );
j = j+1;

% Composition of a CHEBFUN2V with two CHEBFUN3 objects.
f1 = chebfun3(@(x,y,z) x);
f2 = chebfun3(@(x,y,z) y + z);
G = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y, [-1, 1, -2, 2 ]);
H = G(f1, f2);
H_true = chebfun3v(@(x,y,z) x, @(x,y,z) y + z, @(x,y,z) x + y + z);
pass(j) = ( norm(H - H_true) < tol );
j = j+1;

% Composition of a CHEBFUN2V (2 components) with a CHEBFUN2V.
F = chebfun2v(@(x,y) x, @(x,y) y);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y);
H = G(F);
pass(j) = ( norm(H - G) < tol );
j = j+1;

% Composition of a CHEBFUN2V (3 components) with a CHEBFUN2V.
F = chebfun2v(@(x,y) x, @(x,y) y);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y, @(x,y) x.*y);
H = G(F);
pass(j) = ( norm(H - G) < tol );
j = j+1;

% Composition of a CHEBFUN2V (2 components) with a CHEBFUN3V.
F = chebfun3v(@(x,y,z) x, @(x,y,z) y + z);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y, [ -1, 1, -2, 2 ]);
H = G(F);
H_true = chebfun3v(@(x,y,z) x + y + z, @(x,y,z) x - y - z);
pass(j) = ( norm(H - H_true) < tol );
j = j+1;

% Composition of a CHEBFUN2V (3 components) with a CHEBFUN3V.
F = chebfun3v(@(x,y,z) x, @(x,y,z) y + z);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y, @(x,y) x.*y, [ -1, 1, -2, 2 ]);
H = G(F);
H_true = chebfun3v(@(x,y,z) x + y + z, @(x,y,z) x - y - z, @(x,y,z) x.*(y + z));
pass(j) = ( norm(H - H_true) < tol );
j = j+1;

% Composition of a CHEBFUN2V with a CHEBFUN2, which will be interpreted as the
% CHEBFUN2V [real(f); imag(f)].  (Regardless if f is real or complex.)
f = chebfun2(@(x,y) x + 1i*y);
G = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y);
H = G(f);
pass(j) = ( norm(H - G) < tol );
j = j+1;

% With real f:
f = chebfun2(@(x,y) x);
G = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y);
H = G(f);
H_true = chebfun2v(@(x,y) x, @(x,y) 0*x, @(x,y) x);
pass(j) = ( norm(H - H_true) < tol );
j = j+1;

% Composition of a CHEBFUN2V with a CHEBFUN3, which will be interpreted as the
% CHEBFUN3V [real(f); imag(f)].  (Regardless if f is real or complex.)
f = chebfun3(@(x,y,z) x + y + 1i*z);
G = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y, [ -2, 2, -1, 1 ]);
H = G(f);
H_true = chebfun3v(@(x,y,z) x + y, @(x,y,z) z, @(x,y,z) x + y + z);
pass(j) = ( norm(H - H_true) < tol );
j = j+1;

% With real f:
f = chebfun3(@(x,y,z) x + y);
G = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y, [ -2, 2, -2, 2 ]);
H = G(f);
H_true = chebfun3v(@(x,y,z) x + y, @(x,y,z) 0*y, @(x,y,z) x + y);
pass(j) = ( norm(H - H_true) < tol );
j = j+1;

% Test composition with one inf by 2 CHEBFUN:
F = chebfun(@(t) [ t, t ]);
G = chebfun2v(@(x,y) x + y, @(x,y) x);
H = G(F);
t = chebfun(@(t) t);
H_true = [ 2*t, t ];
pass(j) = ( norm(H - H_true) < tol );
j = j+1;

% Test composition with two CHEBFUNs:
f = chebfun(@(t) t);
G = chebfun2v(@(x,y) x + y, @(x,y) x);
H = G(f, f);
H_true = [ 2*f, f ];
pass(j) = ( norm(H - H_true) < tol );
j = j+1;

% Test composition with a complex-valued CHEBFUN:
f = chebfun(@(t) t + 1i*t);
G = chebfun2v(@(x,y) x + y, @(x,y) x);
H = G(f);
t = chebfun(@(t) t);
H_true = [ 2*t, t ];
pass(j) = ( norm(H - H_true) < tol );

end
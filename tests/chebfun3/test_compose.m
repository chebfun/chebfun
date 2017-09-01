function pass = test_compose(pref)
% Check that composition operations are working. 

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) cos(x.*y.*z) + sin(x.*y.*z) + y -.1); 

%% Multiplication
exact = @(x,y,z) (cos(x.*y.*z) + sin(x.*y.*z) + y -.1) .* ...
    sin((x-.1).*(y+.4).*(z+.8));
g = chebfun3(@(x,y,z) exact(x,y,z)); 
[x, y, z] = cheb.xyz;
pass(1) = norm(g - f.*sin((x-.1).*(y+.4).*(z+.8))) < tol; 

%% Sine
exact = @(x,y,z) sin(cos(x.*y.*z) + sin(x.*y.*z) + y -.1); 
g = chebfun3(@(x,y,z) exact(x,y,z));
pass(2) = norm(g - sin(f)) < tol;

g1 = chebfun3(@(x,y,z) g(x,y,z), 'fiberDim', 1); % Construction from a 
% CHEBFUN3 object (as opposed to construction from a function handle)
pass(3) = norm(g1 - sin(f)) < tol;

g2 = chebfun3(@(x,y,z) g(x,y,z), 'fiberDim', 2); 
pass(4) = norm(g2 - sin(f)) < tol;

g3 = chebfun3(@(x,y,z) g(x,y,z), 'fiberDim', 3);
pass(5) = norm(g3 - sin(f)) < tol;


%% Cosine
exact = @(x,y,z) cos(cos(x.*y.*z) + sin(x.*y.*z) + y -.1); 
g = chebfun3(@(x,y,z) exact(x,y,z)); 
pass(6) = norm(g - cos(f)) < tol;

g1 = chebfun3(@(x,y,z) g(x,y,z), 'fiberDim', 1); 
pass(7) = norm(g1 - cos(f)) < tol;

g2 = chebfun3(@(x,y,z) g(x,y,z), 'fiberDim', 2); 
pass(8) = norm(g2 - cos(f)) < tol;

g3 = chebfun3(@(x,y,z) g(x,y,z), 'fiberDim', 3);
pass(9) = norm(g3 - cos(f)) < tol;

%% Sinh
exact = @(x,y,z) sinh(cos(x.*y.*z) + sin(x.*y.*z) + y -.1);
g = chebfun3(@(x,y,z) exact(x,y,z));
pass(10) = norm(g - sinh(f)) < tol;

g1 = chebfun3(@(x,y,z) g(x,y,z), 'fiberDim', 1); 
pass(11) = norm(g1 - sinh(f)) < tol;

g2 = chebfun3(@(x,y,z) g(x,y,z), 'fiberDim', 2); 
pass(12) = norm(g2 - sinh(f)) < tol;

g3 = chebfun3(@(x,y,z) g(x,y,z), 'fiberDim', 3);
pass(13) = norm(g3 - sinh(f)) < tol;


%% Cosh
exact = @(x,y,z) cosh(cos(x.*y.*z) + sin(x.*y.*z) + y -.1);
g = chebfun3(@(x,y,z) exact(x,y,z)); 
pass(14) = norm(g - cosh(f)) < tol;

g1 = chebfun3(@(x,y,z) g(x,y,z), 'fiberDim', 1); 
pass(15) = norm(g1 - cosh(f)) < tol;

g2 = chebfun3(@(x,y,z) g(x,y,z), 'fiberDim', 2); 
pass(16) = norm(g2 - cosh(f)) < tol;

g3 = chebfun3(@(x,y,z) g(x,y,z), 'fiberDim', 3);
pass(17) = norm(g3 - cosh(f)) < tol;

%% Tanh
exact = @(x,y,z) tanh(-(cos(x.*y.*z) + sin(x.*y.*z) + y -.1)); 
g = chebfun3(@(x,y,z) exact(x,y,z)); 
pass(18) = norm(g - tanh(-f)) < tol;

g1 = chebfun3(@(x,y,z) g(x,y,z), 'fiberDim', 1); 
pass(19) = norm(g1 - tanh(-f)) < tol;

g2 = chebfun3(@(x,y,z) g(x,y,z), 'fiberDim', 2); 
pass(20) = norm(g2 - tanh(-f)) < tol;

g3 = chebfun3(@(x,y,z) g(x,y,z), 'fiberDim', 3);
pass(21) = norm(g3 - tanh(-f)) < tol;

%% Exp
exact = @(x,y,z) exp(cos(x.*y.*z) + sin(x.*y.*z) + y -.1); 
g = chebfun3(@(x,y,z) exact(x,y,z)); 
pass(22) = norm(g - exp(f)) < tol;

g1 = chebfun3(@(x,y,z) g(x,y,z), 'fiberDim', 1); 
pass(23) = norm(g1 - exp(f)) < tol;

g2 = chebfun3(@(x,y,z) g(x,y,z), 'fiberDim', 2); 
pass(24) = norm(g2 - exp(f)) < tol;

g3 = chebfun3(@(x,y,z) g(x,y,z), 'fiberDim', 3); 
pass(25) = norm(g3 - exp(f)) < tol;

%% Multiple operations: 
f = chebfun3(@(x,y,z) sin(10*x.*y.*z), [-1 2 -1 1 -3 -1]);
pass(26) = norm(f+f+f-3*f) < 100*tol;

pass(27) = norm(f.*f-f.^2) < tol;

%% Compose a CHEBFUN3 f with a CHEBFUN g (one column):
f = chebfun3(@(x,y,z) x);
g = chebfun(@(t) t.^2);
h = compose(f, g);
pass(28) = ( norm(h - f.^2) < tol );

%% Compose a CHEBFUN3 f with a CHEBFUN g (two columns):
f = chebfun3(@(x,y,z) x);
g = chebfun(@(t) [ t, t.^2 ]);
h = compose(f, g);
pass(29) = ( norm(h - [ f; f.^2 ]) < tol );

%% Compose a complex-valued CHEBFUN3 f with a CHEBFUN2 g:
f = chebfun3(@(x,y,z) x + y + 1i*z);
g = chebfun2(@(x,y) x.^2 + y.^2, [-2, 2, -1, 1]);
h_true = chebfun3(@(x,y,z) (x + y).^2 + z.^2);
h = compose(f, g);
pass(30) = ( norm(h - h_true) < tol );
pass(31) = ~isPeriodicTech(h);

%% Compose a periodic CHEBFUN3 f with a CHEBFUN g
f = chebfun3(@(x,y,z) cos(pi*x), 'trig');
g = chebfun(@(t) t.^2, [-1, 1]);
h = compose(f, g);
pass(32) = isPeriodicTech(h);
G = chebfun(@(t) [ t.^2, cos(t) ], [-1, 1]);
H = compose(f, G);
pass(33) = isPeriodicTech(H);

%% Compose a periodic complex CHEBFUN3 f with a CHEBFUN2 or CHEBFUN2V g:
f = chebfun3(@(x,y,z) exp(1i*pi*x), 'trig');
g = chebfun2(@(x,y) x + y, [-2, 2, -2, 2 ]);
h = compose(f, g);
pass(34) = isPeriodicTech(h);
G = chebfun2v(@(x,y) x, @(x,y) y, [ -2, 2, -2, 2 ]);
H = compose(f, G);
pass(35) = isPeriodicTech(H);

end
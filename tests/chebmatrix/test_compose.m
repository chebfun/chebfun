function pass = test_compose(pref)
% Test COMPOSE().

% Get preferences:
if ( nargin == 0 )
    pref = chebfunpref;
end
tol = 1000*pref.techPrefs.chebfuneps;
j = 1;

% Compose a CHEBMATRIX F with a CHEBFUN2 G:
t = chebfun(@(t) t);
F = [t; t.^2];
g = chebfun2(@(x,y) x + y);
h = compose(F, g); % g(f)
h_true = chebfun(@(t) t + t.^2);
pass(j) = ( norm(h - h_true) < tol );
j = j+1;

% Compose a CHEBMATRIX F with a CHEBFUN2 G and domain not [-1, 1]:
t = chebfun(@(t) t, [0, 2*pi]);
F = [cos(t); sin(t)];
g = chebfun2(@(x,y) x.^2 + y.^2);
h = compose(F, g); % g(f)
h_true = chebfun(@(t) 1 + 0*t, [0, 2*pi]);
pass(j) = ( norm(h - h_true) < tol );
j = j+1;

% Compose a CHEBMATRIX F and a CHEBFUN2V G:
t = chebfun(@(t) t);
F = [t; t.^2];
G = chebfun2v(@(x,y) x + y, @(x,y) y, @(x,y) x);
H = compose(F, G);
H_true = [t + t.^2; t.^2; t];
pass(j) = ( norm(H - H_true) < tol );
j = j+1;

% Compose a CHEBMATRIX F and a CHEBFUN3:
t = chebfun(@(t) t, [0, 2*pi]);
F = [cos(t); sin(t); t];
g = chebfun3(@(x,y,z) x.^2 + y.^2 + z);
h = compose(F, g);
h_true = chebfun(@(t) 1 + t, [0, 2*pi]);
pass(j) = ( norm(h - h_true) < tol );
j = j+1;

% Compose a CHEBMATRIX F and a CHEBFUN3V:
t = chebfun(@(t) t, [0, 2*pi]);
F = [cos(t); sin(t); t];
g = chebfun3v(@(x,y,z) x.^2 + y.^2, @(x,y,z) z);
h = compose(F, g);
h_true = [1 + 0.*t; t];
pass(j) = ( norm(h - h_true) < tol );

end
function pass = test_systemapply
% Test 2x2 systems applied to functions.
% (A Level 3 chebtest.)
% Toby Driscoll, 3 Feb 2014

tol = 1e-14;

d = [-pi,pi];
[Z,I,D] = linop.primitiveOperators(d);
x = chebfun('x',d);

%%
A = [I+2*D^2 -D; D Z];
B = A^2;
u = [ sin(x); exp(x) ];

%%
v = A*u;
w = A*v;
z = (B+3*identity(B))*u;

u1 = u{1}; u2 = u{2};

%%
pass(1) = norm( v{1} - (u1+2*diff(u1,2)-diff(u2))) < 50*tol;
pass(2) = norm( v{2} - (diff(u1))) < 2*tol;
r = w-B*u;
pass = [pass(:); cellfun(@norm,r.blocks) < 1e-10*(tol/eps) ];
r = z - w - 3*u;
pass = [pass(:); cellfun(@norm,r.blocks) < 1e-10*(tol/eps) ];

end

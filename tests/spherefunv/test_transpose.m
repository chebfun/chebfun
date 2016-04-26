function pass = test_transpose( ) 
% Test transpose and ctranspose

tol = 10*chebfunpref().cheb2Prefs.chebfun2eps;

% Test empty spherefunv.
f = spherefunv;
pass(1) = isempty(f');
pass(2) = isempty(f.');

% Test function
f = spherefun(@(x,y,z) cos((x+.1).*y.*z));
% Spherefunv
u = grad(f);
[m,n,p] = size(u');
pass(3) = all( ([isinf(m) isinf(n) p==3]) );
[m,n,p] = size(u.');
pass(4) = all( ([isinf(m) isinf(n) p==3]) );

% Check transpose and ctranspose give the same results for real-valued
% spherefunv objects.
pass(5) = ( norm(u'-u.') < tol );

% Check transpose of transpose if the original u.
v = u';
pass(6) = ( norm(u-v') < tol );

% Check u'*u is a spherefun.
pass(7) = isa(u'*u,'spherefun');

end
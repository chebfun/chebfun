function pass = test_transpose( ) 
% Test transpose and ctranspose

tol = 10*chebfunpref().cheb2Prefs.chebfun2eps;

% Test empty diskfunv.
f = diskfunv;
pass(1) = isempty(f');
pass(2) = isempty(f.');

% Test function
f = diskfun(@(x,y) cos((x+.1).*y));
% diskfunv
u = grad(f);
[m,n,p] = size(u');
pass(3) = all( ([isinf(m) isinf(n) p==2]) );
[m,n,p] = size(u.');
pass(4) = all( ([isinf(m) isinf(n) p==2]) );

% Check transpose and ctranspose give the same results for real-valued
% diskfunv objects.
w = u'-u.';
rng(10); th0 = rand; r0 = rand;
pass(5) = ( norm(w(th0,r0, 'polar')) < tol );

% Check transpose of transpose is the original u.
v = u.';
pass(6) = ( norm(u-v.') < tol );

% Check u'*u is a diskfun.
pass(7) = isa(u'*u,'diskfun');

end
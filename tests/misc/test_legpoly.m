% Test file for LEGPOLY.

function pass = test_legpoly(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

%% Test method 1 on [-1 1]:

% Test 1 confirms that the elements of the array-valued chebfun p 
% are Legendre polynomials: 
p = legpoly(900:1100, [-1 1], 0, 1);
xx = linspace(-1, 1, 10);
P = legendre(900, xx); 
err = norm(feval(p(:,1), xx) - P(1, :), inf);
tol = 5e3*epslevel(p)*vscale(p);
pass(1) = err < tol;


% Test 2 confirms orthogonality:
err = norm(p'*p - diag(diag(p'*p)));
tol = 5e1*epslevel(p)*vscale(p);
pass(2) = err < tol;

% Test 3 confirms othonormality:
p = legpoly(900:1100, [-1 1], 'normalize',1);
err = norm(p'*p - eye(201));
tol = 5e3*epslevel(p)*vscale(p);
pass(3) = err < tol;

%% Test method 1 on [-1 1] (no domain passed):

% Test 4 confirms that the elements of the array-valued chebfun p 
% are Legendre polynomials: 
p = legpoly(900:1100,[-1,1],0,1);
xx = linspace(-1, 1, 10);
P = legendre(900, xx); 
err = norm(feval(p(:,1), xx) - P(1, :), inf);
tol = 5e3*epslevel(p)*vscale(p);
pass(4) = err < tol;

% Test 5 confirms orthogonality:
err = norm(p'*p - diag(diag(p'*p)));
tol = 5e1*epslevel(p)*vscale(p);
pass(5) = err < tol;

% Test 6 confirms othonormality:
p = legpoly(900:1100, 'normalize',1);
err = norm(p'*p - eye(201));
tol = 5e3*epslevel(p)*vscale(p);
pass(6) = err < tol;

%% Test method 1 on [0 10000]:

% Test 7 confirms orthogonality:
p = legpoly(900:1100, [0 10000], 'normalize',1);
err = norm(p'*p - diag(diag(p'*p)));
tol = 1e5*epslevel(p)*vscale(p);
pass(7) = err < tol;

% Test 8 confirms othonormality:
p = legpoly(900:1100, [0 10000], 'normalize',1);
err = norm(p'*p - eye(201));
tol = 1e5*epslevel(p)*vscale(p);
pass(8) = err < tol;

%% Test method 2 on [-1 1]:

% Test 9 confirms that p is a Legendre polynomial: 
p = legpoly(40, [-1 0.2 1]);
xx = linspace(-1, 1, 10);
P = legendre(40, xx);
err = norm(feval(p, xx) - P(1, :), inf);
pass(9) = err < 50*epslevel(p)*vscale(p);

% Test 10 confirms orthogonality:
p = legpoly(1:100, [-1 -0.2 0.3 1]);
err = norm(p'*p - diag(diag(p'*p)));
pass(10) = err < 10*epslevel(p)*vscale(p);

% Test 11 confirms othonormality:
p = legpoly(1:100, [-1 0.145 1], 'normalize');
err = norm(p'*p - eye(100));
tol = 100*epslevel(p)*vscale(p);
pass(11) = err < tol;


%% Test method 2 on [-1 1] (no domain passed):

% Test 12 confirms that p is a Legendre polynomial: 
p = legpoly(40);
xx = linspace(-1, 1, 10);
P = legendre(40, xx);
err = norm(feval(p, xx) - P(1, :), inf);
tol = 50*epslevel(p)*vscale(p);
pass(12) = err < tol;

% Test 13 confirms orthogonality:
p = legpoly(1:100);
err = norm(p'*p - diag(diag(p'*p)));
tol = 10*epslevel(p)*vscale(p);
pass(13) = err < tol;

% Test 14 confirms othonormality:
p = legpoly(1:100, 'normalize');
err = norm(p'*p - eye(100));
tol = 100*epslevel(p)*vscale(p);
pass(14) = err < tol;

%% Test method 2 on [0 10000]:

% Test 15 confirms orthogonality:
p = legpoly(1:100, [0 155 3333 10000]);
err = norm(p'*p - diag(diag(p'*p)));
pass(15) = err < 5e4*epslevel(p)*vscale(p);

% Test 16 confirms othonormality:
p = legpoly(1:100, [0 3333 10000], 'normalize');
err = norm(p'*p - eye(100));
tol = 1e4*epslevel(p)*vscale(p);
pass(16) = err < tol;

%% Test method 3 on [-1 1]:

% Test 17 confirms that p is a Legendre polynomial: 
p = legpoly(1500, [-1 -0.5 -0.3 1]);
xx = linspace(-1, 1, 10);
P = legendre(1500, xx);
err = norm(feval(p, xx) - P(1, :), inf);
tol = 1e4*epslevel(p)*vscale(p);
pass(17) = err < tol;

% Test 18 confirms normaliztion:
p = legpoly(1500, [-1 0 1], 'normalize');
err = norm(p'*p-1);
pass(18) = err < 1e4*epslevel(p)*vscale(p);

%% Test method 3 on [-1 1] (no domain passed):

% Test 19 confirms that p is a Legendre polynomial: 
p = legpoly(1500);
xx = linspace(-1, 1, 10);
P = legendre(1500, xx);
err = norm(feval(p, xx) - P(1, :), inf);
tol = 1e4*epslevel(p)*vscale(p);
pass(19) = err < tol;

% Test 20 confirms normaliztion:
p = legpoly(1500, 'normalize');
err = norm(p'*p-1);
tol = 1e4*epslevel(p)*vscale(p);
pass(20) = err < tol;

% Test 21 & 22 confirms vectorized version of METHOD 3:
p = legpoly([1000 1500]);
xx = linspace(-1, 1, 10);
P1 = legendre(1000, xx);
P2 = legendre(1500, xx);
err1 = norm(feval(p(:,1), xx) - P1(1, :), inf);
err2 = norm(feval(p(:,2), xx) - P2(1, :), inf);
tol = 1e4*epslevel(p)*vscale(p);
pass(21) = err1 < tol;
pass(22) = err2 < tol;

% Test 23 confirms normaliztion for vectorized version of METHOD 3:
p = legpoly([1000 1500], 'normalize');
err = norm(p'*p-eye(2));
tol = 1e4*epslevel(p)*vscale(p);
pass(23) = err < tol;

%% Test method 2 on [-1 1] with a row input N:

% Test 24 confirms orthogonality:
p = legpoly([1;2;3;4;5;6;7;8;9;10], [-1 0.2 1]);
err = norm(p*p' - diag(diag(p*p')));
tol = 5*epslevel(p)*vscale(p);
pass(24) = err < tol;

% Test 25 confirms othonormality:
p = legpoly([1;2;3;4;5;6;7;8;9;10], [-1 0.6 1], 'normalize');
err = norm(p*p' - eye(10));
tol = 10*epslevel(p)*vscale(p);
pass(25) = err < tol;

% Test 26 checks empty case:
p = legpoly([], [-1 0.66 1], 'normalize');
pass(26) = isempty(p);

end
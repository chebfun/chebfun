% Test file for LEGPOLY.

function pass = test_legpoly(pref)

if ( nargin == 0 )
    pref = chebpref();
end

%% Test method 1 on [-1 1]:

% Test 1 confirms that the elements of the array-valued chebfun p 
% are Legendre polynomials: 
p = legpoly(900:1100, [-1 0 1]);
xx = linspace(-1, 1, 1000);
P = legendre(900, xx); 
err = norm(feval(p(:,1), xx) - P(1, :), inf);
pass(1) = err < 5e3*epslevel(p)*vscale(p);

% Test 2 confirms orthogonality:
err = norm(p'*p - diag(diag(p'*p)));
pass(2) = err < 5e1*epslevel(p)*vscale(p);

% Test 3 confirms othonormality:
p = legpoly(900:1100, [-1 0 1], 'normalize');
err = norm(p'*p - eye(201));
pass(3) = err < 5e3*epslevel(p)*vscale(p);

%% Test method 1 on [-1 1] (no domain passed):

% Test 4 confirms that the elements of the array-valued chebfun p 
% are Legendre polynomials: 
p = legpoly(900:1100);
xx = linspace(-1, 1, 1000);
P = legendre(900, xx); 
err = norm(feval(p(:,1), xx) - P(1, :), inf);
pass(4) = err < 5e3*epslevel(p)*vscale(p);

% Test 5 confirms orthogonality:
err = norm(p'*p - diag(diag(p'*p)));
pass(5) = err < 5e1*epslevel(p)*vscale(p);

% Test 6 confirms othonormality:
p = legpoly(900:1100, 'normalize');
err = norm(p'*p - eye(201));
pass(6) = err < 5e3*epslevel(p)*vscale(p);

%% Test method 1 on [0 10000]:

% Test 7 confirms orthogonality:
p = legpoly(900:1100, [0 100 10000], 'normalize');
err = norm(p'*p - diag(diag(p'*p)));
pass(7) = err < 5e4*epslevel(p)*vscale(p);

% Test 8 confirms othonormality:
p = legpoly(900:1100, [0 50 10000], 'normalize');
err = norm(p'*p - eye(201));
pass(8) = err < 5e4*epslevel(p)*vscale(p);

%% Test method 2 on [-1 1]:

% Test 9 confirms that p is a Legendre polynomial: 
p = legpoly(40, [-1 0.2 1]);
xx = linspace(-1, 1, 1000);
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
pass(11) = err < 100*epslevel(p)*vscale(p);


%% Test method 2 on [-1 1] (no domain passed):

% Test 12 confirms that p is a Legendre polynomial: 
p = legpoly(40);
xx = linspace(-1, 1, 1000);
P = legendre(40, xx);
err = norm(feval(p, xx) - P(1, :), inf);
pass(12) = err < 50*epslevel(p)*vscale(p);

% Test 13 confirms orthogonality:
p = legpoly(1:100);
err = norm(p'*p - diag(diag(p'*p)));
pass(13) = err < 10*epslevel(p)*vscale(p);

% Test 14 confirms othonormality:
p = legpoly(1:100, 'normalize');
err = norm(p'*p - eye(100));
pass(14) = err < 100*epslevel(p)*vscale(p);

%% Test method 2 on [0 10000]:

% Test 15 confirms orthogonality:
p = legpoly(1:100, [0 155 3333 10000]);
err = norm(p'*p - diag(diag(p'*p)));
pass(15) = err < 5e4*epslevel(p)*vscale(p);

% Test 16 confirms othonormality:
p = legpoly(1:100, [0 3333 10000], 'normalize');
err = norm(p'*p - eye(100));
pass(16) = err < 5e3*epslevel(p)*vscale(p);

%% Test method 3 on [-1 1]:

% Test 17 confirms that p is a Legendre polynomial: 
p = legpoly(1500, [-1 -0.5 -0.3 1]);
xx = linspace(-1, 1, 1000);
P = legendre(1500, xx);
err = norm(feval(p, xx) - P(1, :), inf);
pass(17) = err < 1e4*epslevel(p)*vscale(p);

% Test 18 confirms normaliztion:
p = legpoly(1500, [-1 0 1], 'normalize');
err = norm(p'*p-1);
pass(18) = err < 1e4*epslevel(p)*vscale(p);

%% Test method 3 on [-1 1] (no domain passed):

% Test 19 confirms that p is a Legendre polynomial: 
p = legpoly(1500);
xx = linspace(-1, 1, 1000);
P = legendre(1500, xx);
err = norm(feval(p, xx) - P(1, :), inf);
pass(19) = err < 1e4*epslevel(p)*vscale(p);

% Test 20 confirms normaliztion:
p = legpoly(1500, 'normalize');
err = norm(p'*p-1);
pass(20) = err < 1e4*epslevel(p)*vscale(p);

%% Test method 2 on [-1 1] with a row input N:

% Test 21 confirms orthogonality:
p = legpoly([1;2;3;4;5;6;7;8;9;10], [-1 0.2 1]);
err = norm(p*p' - diag(diag(p*p')));
pass(21) = err < 5*epslevel(p)*vscale(p);

% Test 22 confirms othonormality:
p = legpoly([1;2;3;4;5;6;7;8;9;10], [-1 0.6 1], 'normalize');
err = norm(p*p' - eye(10));
pass(22) = err < 10*epslevel(p)*vscale(p);

% Test 23 checks empty case:
p = legpoly([], [-1 0.66 1], 'normalize');
pass(23) = isempty(p);

end
function pass = test_diffmat(pref)

if ( nargin == 0 )
    pref = cheboppref();
end

tol = 1e-10;

c1 = colloc1();
c2 = colloc2();

N = 5;
Da = diffmat(N);
Db = colloc2.diffmat(N);
err = norm(Da - Db);
pass(1) = err < tol;

N = 5;
Da = diffmat(N, [-2 2]);
Db = c2.diffmat(N)/2;
err = norm(Da - Db);
pass(2) = err < tol;

N = 5;
Da = diffmat(N, 2, [-2 2]);
Db = c2.diffmat(N,2)/4;
err = norm(Da - Db);
pass(3) = err < tol;

N = 5;
Da = diffmat(N, c1);
Db = c1.diffmat(N,1);
err = norm(Da - Db);
pass(4) = err < tol;

N = 5;
Da = diffmat(N, 3, [-.5 .5], 'colloc1');
Db = c1.diffmat(N,3)*8;
err = norm(Da - Db);
pass(5) = err < tol;

p = cheboppref();
cheboppref.setDefaults('discretization', @colloc1);
try
    N = 5;
    Da = diffmat(N, 3);
    Db = c2.diffmat(N,3);
    err = norm(Da - Db);
    pass(6) = err < tol;
catch ME
    chebopref.setDefaults(p);
    rethrow(ME)
end
cheboppref.setDefaults(p);

%% periodic function:
N = 6;
dom = [0 2*pi];
op = @(x)sin(x);
opp = @(x)-cos(x);
x = fourtech.fourpts(N);
x = diff(dom)/2*(x-dom(1));
f = op(x);
fp_exact = opp(x);
D = diffmat(N, 3, 'periodic', dom);
fp = D*f;
err = norm(fp-fp_exact, inf);
pass(6) = ( err < tol );

%% Rectangular differentiation matrices:

%% 1st-kind grid -> 1st-kind grid:

% P = N-M & [-1 1];
N = 5;
p = 2;
x = chebpts(N, 1);
y = chebpts(N-p, 1);
op = @(x)x.^3 + 1;
opp = @(x)6*x;
f = op(x);
fp_exact = opp(y);
D = diffmat([N-p N], p, [-1 1], 'colloc1');
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(7) = ( err < tol );

% P = N-M & arbitrary domain:
dom = [-2 7];
x = chebpts(N, dom, 1);
y = chebpts(N-p, dom, 1);
f = op(x);
fp_exact = opp(y);
D = diffmat([N-p N], p, dom, 'colloc1');
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(8) = ( err < tol );

% Other syntax:
DD = diffmat([N -p], p, dom, 'colloc1');
err = norm(D-DD, inf);
pass(9) = ( err < tol );

% Other syntax 2:
DD = diffmat(N, p, 'rect', dom, 'colloc1');
err = norm(D-DD, inf);
pass(10) = ( err < tol );

% Other syntax 3:
DD = diffmat(N, p, dom, 'rect', 'colloc1');
err = norm(D-DD, inf);
pass(11) = ( err < tol );

% P ~= N-M & arbitrary domain:
dom = [-2 7];
M = N-p-1;
x = chebpts(N, dom, 1);
y = chebpts(M, dom, 1);
f = op(x);
fp_exact = opp(y);
D = diffmat([M N], p, dom, 'colloc1');
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(12) = ( err < tol );

%% 2nd-kind grid -> 1st-kind grid:

% P = N-M & [-1 1];
N = 5;
p = 2;
x = chebpts(N);
y = chebpts(N-p, 1);
f = op(x);
fp_exact = opp(y);
D = diffmat([N-p N], p);
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(13) = ( err < tol );

% P = N-M & arbitrary domain:
dom = [-2 7];
x = chebpts(N, dom);
y = chebpts(N-p, dom, 1);
f = op(x);
fp_exact = opp(y);
D = diffmat([N-p N], p, dom);
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(14) = ( err < tol );

% Other syntax:
DD = diffmat([N -p], p, dom);
err = norm(D-DD, inf);
pass(15) = ( err < tol );

% Other syntax 2:
DD = diffmat(N, p, 'rect', dom);
err = norm(D-DD, inf);
pass(16) = ( err < tol );

% Other syntax 3:
DD = diffmat(N, p, dom, 'rect');
err = norm(D-DD, inf);
pass(17) = ( err < tol );

% P ~= N-M & arbitrary domain:
dom = [-2 7];
M = N-p-1;
x = chebpts(N, dom);
y = chebpts(M, dom, 1);
f = op(x);
fp_exact = opp(y);
D = diffmat([M N], p, dom);
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(18) = ( err < tol );

%% Boundary conditions:

%% 1st-kind grid -> 1st-kind grid:

% 1st-order problem: u' = exp(x); u(-2) = exp(-2):
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 1;
M = N - p;
x = chebpts(N, dom, 1);
y = chebpts(M, dom, 1);
f = op(x);
fp = op(y);
rhs = [op(dom(1)); fp];
D = diffmat([M N], p, dom, 'colloc1', {'d'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(19) = ( err < tol );

% 1st-order problem: u' = exp(x); u(7) = exp(7):
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 1;
M = N - p;
x = chebpts(N, dom, 1);
y = chebpts(M, dom, 1);
f = op(x);
fp = op(y);
rhs = [fp; op(dom(2))];
D = diffmat([M N], p, dom, 'colloc1', {}, {'d'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(20) = ( err < tol );


% 2nd-order problem: u" = exp(x); u(-2) = exp(-2); u(7) = exp(7);
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 2;
M = N - p;
x = chebpts(N, dom, 1);
y = chebpts(M, dom, 1);
f = op(x);
fp = op(y);
rhs = [op(dom(1)); fp; op(dom(2))];
D = diffmat([M N], p, dom, 'colloc1', {'d'}, {'d'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(21) = ( err < 1e1*tol );

% 2nd-order problem: u" = exp(x); u'(2) = exp(2); u(2) = exp(2);
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 2;
M = N - p;
x = chebpts(N, dom, 1);
y = chebpts(M, dom, 1);
f = op(x);
fp = op(y);
rhs = [op(dom(1)); op(dom(1)); fp];
D = diffmat([M N], p, dom, 'colloc1', {'d' 'n'}, {});
ff = D\rhs;
err = norm(ff-f, inf);
pass(22) = ( err < 1e2*tol );

% 2nd-order problem: u" = exp(x); u'(7) = exp(7); u(7) = exp(7);
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 2;
M = N - p;
x = chebpts(N, dom, 1);
y = chebpts(M, dom, 1);
f = op(x);
fp = op(y);
rhs = [fp; op(dom(2)); op(dom(2))];
D = diffmat([M N], p, dom, 'colloc1', {}, {'n' 'd'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(23) = ( err < 1e1*tol );

% 2nd-order problem: u" = exp(x); u(-2) = exp(-2); u'(7) = exp(7);
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 2;
M = N - p;
x = chebpts(N, dom, 1);
y = chebpts(M, dom, 1);
f = op(x);
fp = op(y);
rhs = [op(dom(1)); fp; op(dom(2))];
D = diffmat([M N], p, dom, 'colloc1', {'d'}, {'n'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(24) = ( err < 1e1*tol );

% 2nd-order problem: u" = exp(x); u'(-2) = exp(-2); u(7) = exp(7);
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 2;
M = N - p;
x = chebpts(N, dom, 1);
y = chebpts(M, dom, 1);
f = op(x);
fp = op(y);
rhs = [op(dom(1)); fp; op(dom(2))];
D = diffmat([M N], p, dom, 'colloc1', {'n'}, {'d'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(25) = ( err < 1e2*tol );

% 2nd-order problem: u" = exp(x); sum(u) = exp(7)-exp(-2); u'(-2) = exp(-2);
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 2;
M = N - p;
x = chebpts(N, dom, 1);
y = chebpts(M, dom, 1);
f = op(x);
fp = op(y);
rhs = [op(dom(1)); fp; op(dom(2))-op(dom(1))];
D = diffmat([M N], p, dom, 'colloc1', {'n'}, {'s'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(26) = ( err < 1e1*tol );

%% 2nd-kind grid -> 1st-kind grid:

% 1st-order problem: u' = exp(x); u(-2) = exp(-2):
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 1;
M = N - p;
x = chebpts(N, dom);
y = chebpts(M, dom, 1);
f = op(x);
fp = op(y);
rhs = [op(dom(1)); fp];
D = diffmat([N -p], p, dom, 'colloc2', {'d'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(27) = ( err < tol );

% 1st-order problem: u' = exp(x); u(7) = exp(7):
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 1;
M = N - p;
x = chebpts(N, dom);
y = chebpts(M, dom, 1);
f = op(x);
fp = op(y);
rhs = [fp; op(dom(2))];
D = diffmat([N -p], p, dom, 'colloc2', {}, {'d'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(28) = ( err < tol );

% 2nd-order problem: u" = exp(x); u(-2) = exp(-2); u(7) = exp(7);
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 2;
M = N - p;
x = chebpts(N, dom);
y = chebpts(M, dom, 1);
f = op(x);
fp = op(y);
rhs = [op(dom(1)); fp; op(dom(2))];
D = diffmat([M N], p, dom, 'colloc2', {'d'}, {'d'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(29) = ( err < 1e1*tol );

% 2nd-order problem: u" = exp(x); u'(2) = exp(2); u(2) = exp(2);
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 2;
M = N - p;
x = chebpts(N, dom);
y = chebpts(M, dom, 1);
f = op(x);
fp = op(y);
rhs = [op(dom(1)); op(dom(1)); fp];
D = diffmat([M N], p, dom, 'colloc2', {'d' 'n'}, {});
ff = D\rhs;
err = norm(ff-f, inf);
pass(30) = ( err < 1e2*tol );

% 2nd-order problem: u" = exp(x); u'(7) = exp(7); u(7) = exp(7);
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 2;
M = N - p;
x = chebpts(N, dom);
y = chebpts(M, dom, 1);
f = op(x);
fp = op(y);
rhs = [fp; op(dom(2)); op(dom(2))];
D = diffmat([M N], p, dom, 'colloc2', {}, {'n' 'd'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(31) = ( err < 1e2*tol );

% 2nd-order problem: u" = exp(x); u(-2) = exp(-2); u'(7) = exp(7);
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 2;
M = N - p;
x = chebpts(N, dom);
y = chebpts(M, dom, 1);
f = op(x);
fp = op(y);
rhs = [op(dom(1)); fp; op(dom(2))];
D = diffmat([M N], p, dom, 'colloc2', {'d'}, {'n'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(32) = ( err < 1e2*tol );

% 2nd-order problem: u" = exp(x); u'(-2) = exp(-2); u(7) = exp(7);
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 2;
M = N - p;
x = chebpts(N, dom);
y = chebpts(M, dom, 1);
f = op(x);
fp = op(y);
rhs = [op(dom(1)); fp; op(dom(2))];
D = diffmat([M N], p, dom, 'colloc2', {'n'}, {'d'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(33) = ( err < 1e2*tol );

% 2nd-order problem: u" = exp(x); sum(u) = exp(7)-exp(-2); u(7) = exp(7);
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 2;
M = N - p;
x = chebpts(N, dom);
y = chebpts(M, dom, 1);
f = op(x);
fp = op(y);
rhs = [fp; op(dom(2))-op(dom(1)); op(dom(2))];
D = diffmat([M N], p, dom, 'colloc2', {}, {'s' 'd'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(34) = ( err < 1e1*tol );

end
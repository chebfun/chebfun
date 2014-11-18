function pass = test_diffmat(pref)

if ( nargin == 0 )
    pref = cheboppref();
end

tol = 1e-10;

c1 = chebcolloc1();
c2 = chebcolloc2();

N = 5;
Da = diffmat(N);
Db = c2.diffmat(N);
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
Da = diffmat(N, 'chebkind1');
Db = c1.diffmat(N, 1);
err = norm(Da - Db);
pass(4) = err < tol;

N = 5;
Da = diffmat(N, 3, [-.5 .5], 'chebkind1');
Db = c1.diffmat(N, 3)*8;
err = norm(Da - Db);
pass(5) = err < tol;

p = cheboppref();
cheboppref.setDefaults('discretization', @chebcolloc1);
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

%% Legendre square case:
N = 6;
p = 2;
x = legpts(N);
op = @(x)x.^3 + 1;
opp = @(x)6*x;
f = op(x);
fp_exact = opp(x);
D = diffmat(N, p, 'leg');
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(7) = ( err < tol );

%% periodic function:
N = 6;
dom = [0 2*pi];
op = @(x)sin(x);
opp = @(x)-cos(x);
x = trigtech.trigpts(N);
x = diff(dom)/2*(x-dom(1));
f = op(x);
fp_exact = opp(x);
D = diffmat(N, 3, 'periodic', dom);
fp = D*f;
err = norm(fp-fp_exact, inf);
pass(8) = ( err < tol );

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
D = diffmat([N-p N], p, [-1 1], 'chebkind1');
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(9) = ( err < tol );

% P = N-M & arbitrary domain:
dom = [-2 7];
x = chebpts(N, dom, 1);
y = chebpts(N-p, dom, 1);
f = op(x);
fp_exact = opp(y);
D = diffmat([N-p N], p, dom, 'chebkind1');
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(10) = ( err < tol );

% Other syntax:
DD = diffmat([N -p], p, dom, 'chebkind1');
err = norm(D-DD, inf);
pass(11) = ( err < tol );

% Other syntax 2:
DD = diffmat(N, p, 'rect', dom, 'chebkind1');
err = norm(D-DD, inf);
pass(12) = ( err < tol );

% Other syntax 3:
DD = diffmat(N, p, dom, 'rect', 'chebkind1');
err = norm(D-DD, inf);
pass(13) = ( err < tol );

% P ~= N-M & arbitrary domain:
dom = [-2 7];
M = N-p-1;
x = chebpts(N, dom, 1);
y = chebpts(M, dom, 1);
f = op(x);
fp_exact = opp(y);
D = diffmat([M N], p, dom, 'chebkind1');
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(14) = ( err < tol );

%% 2nd-kind grid -> 1st-kind grid:

% P = N-M & [-1 1];
N = 5;
p = 2;
x = chebpts(N);
y = chebpts(N-p, 1);
f = op(x);
fp_exact = opp(y);
D = diffmat([N-p N], p, 'chebkind2', 'chebkind1');
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(15) = ( err < tol );

% P = N-M & arbitrary domain:
dom = [-2 7];
x = chebpts(N, dom);
y = chebpts(N-p, dom, 1);
f = op(x);
fp_exact = opp(y);
D = diffmat([N-p N], p, dom, 'chebkind2', 'chebkind1');
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(16) = ( err < tol );

% Other syntax:
DD = diffmat([N -p], p, dom, 'chebkind2', 'chebkind1');
err = norm(D-DD, inf);
pass(17) = ( err < tol );

% Other syntax 2:
DD = diffmat(N, p, 'rect', dom, 'chebkind2', 'chebkind1');
err = norm(D-DD, inf);
pass(18) = ( err < tol );

% Other syntax 3:
DD = diffmat(N, p, dom, 'rect', 'chebkind2', 'chebkind1');
err = norm(D-DD, inf);
pass(19) = ( err < tol );

% P ~= N-M & arbitrary domain:
dom = [-2 7];
M = N-p-1;
x = chebpts(N, dom);
y = chebpts(M, dom, 1);
f = op(x);
fp_exact = opp(y);
D = diffmat([M N], p, dom, 'chebkind2', 'chebkind1');
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(20) = ( err < tol );

%% Legendre grid -> Legendre grid:

% P = N-M & [-1 1];
N = 6;
p = 2;
x = legpts(N);
y = legpts(N-p);
f = op(x);
fp_exact = opp(y);
D = diffmat([N-p N], p, 'leg', 'leg');
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(21) = ( err < tol );

% P = N-M & arbitrary domain:
dom = [-2 7];
x = legpts(N, dom);
y = legpts(N-p, dom);
f = op(x);
fp_exact = opp(y);
D = diffmat([N-p N], p, dom, 'leg', 'leg');
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(22) = ( err < tol );

% Other syntax:
DD = diffmat([N -p], p, dom, 'leg');
err = norm(D-DD, inf);
pass(23) = ( err < tol );

% Other syntax 2:
DD = diffmat(N, p, 'rect', dom, 'leg');
err = norm(D-DD, inf);
pass(24) = ( err < tol );

% Other syntax 3:
DD = diffmat(N, p, dom, 'rect', 'leg');
err = norm(D-DD, inf);
pass(25) = ( err < tol );

% P ~= N-M & arbitrary domain:
dom = [-2 7];
M = N-p-1;
x = legpts(N, dom);
y = legpts(M, dom);
f = op(x);
fp_exact = opp(y);
D = diffmat([M N], p, dom, 'leg');
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(26) = ( err < tol );

%% Legendre grid -> 1st-kind grid:

% P = N-M & arbitrary domain:
dom = [-2 7];
x = legpts(N, dom);
y = chebpts(N-p, dom, 1);
f = op(x);
fp_exact = opp(y);
D = diffmat([N-p N], p, dom, 'leg', 'chebkind1');
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(27) = ( err < tol );

%% 1st-kind grid -> Legendre grid:

% P = N-M & arbitrary domain:
dom = [-2 7];
x = chebpts(N, dom, 1);
y = legpts(N-p, dom);
f = op(x);
fp_exact = opp(y);
D = diffmat([N-p N], p, dom, 'chebkind1', 'leg');
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(28) = ( err < tol );

%% 2nd-kind grid -> Legendre grid:

% P = N-M & arbitrary domain:
dom = [-2 7];
x = chebpts(N, dom);
y = legpts(N-p, dom);
f = op(x);
fp_exact = opp(y);
D = diffmat([N-p N], p, dom, 'chebkind2', 'leg');
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(29) = ( err < tol );

%% 1st-kind grid -> 2nd-kind grid:

% P = N-M & arbitrary domain:
dom = [-2 7];
x = chebpts(N, dom, 1);
y = chebpts(N-p, dom);
f = op(x);
fp_exact = opp(y);
D = diffmat([N-p N], p, dom, 'chebkind1', 'chebkind2');
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(30) = ( err < tol );

%% 2nd-kind grid -> 2nd-kind grid:

% P = N-M & arbitrary domain:
dom = [-2 7];
x = chebpts(N, dom);
y = chebpts(N-p, dom);
f = op(x);
fp_exact = opp(y);
D = diffmat([N-p N], p, dom, 'chebkind2', 'chebkind2');
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(31) = ( err < tol );

%% Legendre grid -> 2nd-kind grid:

% P = N-M & arbitrary domain:
dom = [-2 7];
x = legpts(N, dom);
y = chebpts(N-p, dom);
f = op(x);
fp_exact = opp(y);
D = diffmat([N-p N], p, dom, 'leg', 'chebkind2');
fp = D*f;
err = norm(fp_exact - fp, inf);
pass(32) = ( err < tol );

%%%%%%%%%%%%%%%%%%%%% Boundary conditions: %%%%%%%%%%%%%%%%%%%%%%%%%%%

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
D = diffmat([M N], p, dom, 'chebkind1', 'dirichlet');
ff = D\rhs;
err = norm(ff-f, inf);
pass(33) = ( err < tol );

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
D = diffmat([M N], p, dom, 'chebkind1', [], 'dirichlet');
ff = D\rhs;
err = norm(ff-f, inf);
pass(34) = ( err < tol );

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
D = diffmat([M N], p, dom, 'chebkind1', {'dirichlet'}, {'dirichlet'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(35) = ( err < 1e1*tol );

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
D = diffmat([M N], p, dom, 'chebkind1', {'dirichlet' 'neumann'}, {});
ff = D\rhs;
err = norm(ff-f, inf);
pass(36) = ( err < 1e2*tol );

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
D = diffmat([M N], p, dom, 'chebkind1', {}, {'neumann' 'dirichlet'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(37) = ( err < 5e1*tol );

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
D = diffmat([M N], p, dom, 'chebkind1', {'dirichlet'}, {'neumann'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(38) = ( err < 2e1*tol );

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
D = diffmat([M N], p, dom, 'chebkind1', {'neumann'}, {'dirichlet'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(39) = ( err < 1e2*tol );

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
D = diffmat([M N], p, dom, 'chebkind1', {'neumann'}, {'sum'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(40) = ( err < 5e1*tol );

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
D = diffmat([N -p], p, dom, 'chebkind2', 'chebkind1', {'dirichlet'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(41) = ( err < tol );

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
D = diffmat([N -p], p, dom, 'chebkind2', 'chebkind1', {}, {'dirichlet'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(42) = ( err < tol );

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
D = diffmat([M N], p, dom, 'chebkind2', 'chebkind1', 'dirichlet', {'dirichlet'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(43) = ( err < 2e1*tol );

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
D = diffmat([M N], p, dom, 'chebkind2', 'chebkind1', {'dirichlet' 'neumann'}, {});
ff = D\rhs;
err = norm(ff-f, inf);
pass(44) = ( err < 1e2*tol );

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
D = diffmat([M N], p, dom, 'chebkind2', 'chebkind1', [], {'neumann' 'dirichlet'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(45) = ( err < 1e2*tol );

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
D = diffmat([M N], p, dom, 'chebkind2', 'chebkind1', {'dirichlet'}, {'neumann'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(46) = ( err < 1e2*tol );

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
D = diffmat([M N], p, dom, 'chebkind2', 'chebkind1', {'neumann'}, {'dirichlet'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(47) = ( err < 1e3*tol );

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
D = diffmat([M N], p, dom, 'chebkind2', 'chebkind1', {}, {'sum' 'dirichlet'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(48) = ( err < 1e1*tol );

%% legendre grid -> legendre grid:

% 1st-order problem: u' = exp(x); u(-2) = exp(-2):
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 1;
M = N - p;
x = legpts(N, dom);
y = legpts(M, dom);
f = op(x);
fp = op(y);
rhs = [op(dom(1)); fp];
D = diffmat([N -p], p, dom, 'leg', 'dirichlet');
ff = D\rhs;
err = norm(ff-f, inf);
pass(49) = ( err < tol );

% 1st-order problem: u' = exp(x); u(7) = exp(7):
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 1;
M = N - p;
x = legpts(N, dom);
y = legpts(M, dom);
f = op(x);
fp = op(y);
rhs = [fp; op(dom(2))];
D = diffmat([N -p], p, dom, 'leg', [], 'dirichlet');
ff = D\rhs;
err = norm(ff-f, inf);
pass(50) = ( err < tol );

% 2nd-order problem: u" = exp(x); u(-2) = exp(-2); u(7) = exp(7);
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 2;
M = N - p;
x = legpts(N, dom);
y = legpts(M, dom);
f = op(x);
fp = op(y);
rhs = [op(dom(1)); fp; op(dom(2))];
D = diffmat([M N], p, dom, 'leg', 'dirichlet', 'dirichlet');
ff = D\rhs;
err = norm(ff-f, inf);
pass(51) = ( err < 1e1*tol );

% 2nd-order problem: u" = exp(x); u'(2) = exp(2); u(2) = exp(2);
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 2;
M = N - p;
x = legpts(N, dom);
y = legpts(M, dom);
f = op(x);
fp = op(y);
rhs = [op(dom(1)); op(dom(1)); fp];
D = diffmat([M N], p, dom, 'leg', {'dirichlet' 'neumann'}, {});
ff = D\rhs;
err = norm(ff-f, inf);
pass(52) = ( err < 1e2*tol );

% 2nd-order problem: u" = exp(x); u'(7) = exp(7); u(7) = exp(7);
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 2;
M = N - p;
x = legpts(N, dom);
y = legpts(M, dom);
f = op(x);
fp = op(y);
rhs = [fp; op(dom(2)); op(dom(2))];
D = diffmat([M N], p, dom, 'leg', [], {'neumann' 'dirichlet'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(53) = ( err < 1e2*tol );

% 2nd-order problem: u" = exp(x); u(-2) = exp(-2); u'(7) = exp(7);
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 2;
M = N - p;
x = legpts(N, dom);
y = legpts(M, dom);
f = op(x);
fp = op(y);
rhs = [op(dom(1)); fp; op(dom(2))];
D = diffmat([M N], p, dom, 'leg', 'dirichlet', 'neumann');
ff = D\rhs;
err = norm(ff-f, inf);
pass(54) = ( err < 1e2*tol );

% 2nd-order problem: u" = exp(x); u'(-2) = exp(-2); u(7) = exp(7);
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 2;
M = N - p;
x = legpts(N, dom);
y = legpts(M, dom);
f = op(x);
fp = op(y);
rhs = [op(dom(1)); fp; op(dom(2))];
D = diffmat([M N], p, dom, 'leg', {'neumann'}, {'dirichlet'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(55) = ( err < 1e2*tol );

% 2nd-order problem: u" = exp(x); sum(u) = exp(7)-exp(-2); u(7) = exp(7);
op = @(x)exp(x);
dom = [-2 7];
N = 30;
p = 2;
M = N - p;
x = legpts(N, dom);
y = legpts(M, dom);
f = op(x);
fp = op(y);
rhs = [fp; op(dom(2))-op(dom(1)); op(dom(2))];
D = diffmat([M N], p, dom, 'leg', [], {'sum' 'dirichlet'});
ff = D\rhs;
err = norm(ff-f, inf);
pass(56) = ( err < 1e1*tol );

%% Test on symmetry:
D = diffmat(3); 
d = D + D(end:-1:1,end:-1:1);
pass(57) = ~norm(d, inf);
D = diffmat(3, 'chebkind1'); 
d = D + D(end:-1:1,end:-1:1);
pass(58) = ~norm(d, inf);
D = diffmat(3, 'leg'); 
d = D + D(end:-1:1,end:-1:1);
pass(59) = ~norm(d, inf);
D = diffmat(4); 
d = D + D(end:-1:1,end:-1:1);
pass(60) = ~norm(d, inf);
D = diffmat(4, 'chebkind1'); 
d = D + D(end:-1:1,end:-1:1);
pass(61) = ~norm(d, inf);
D = diffmat(4, 'leg'); 
d = D + D(end:-1:1,end:-1:1);
pass(62) = ~norm(d, inf);

end

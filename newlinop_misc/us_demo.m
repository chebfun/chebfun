ccc
format short

%% Differentiation matrix

%%
% Define a differential operator on [-1, 1]:
dom = [-1, 1];
D = linop.diff(dom)

%%
% Return its 5x5 US discretization:
n = 9;
M = discretize(D, n, dom, @blockUS)
full(M)

%%
% It takes this form because $$\frac{d}{dx}T_k(x) = kU_{k-1}(x).$$

%%
% M*a = b maps a function $u(x) = \sum a_k T_k(x)$ to $u'(x) = \sum b_k U_k(x)$.
f = chebfun(@cos, n);
a = flipud(f.coeffs);
b = M*a

%%
% We then need the conversion matrix S so that $$u'(x) = \sum b_k U_k(x) = \sum
% c_k T_k(x)$$, where S*c = b.

S = blockUS.convertmat(n, 0, 0)
full(S)
c = S\b

%%
% We can compare this to computing the deriviative in CHEBFUN:
fp = diff(f);
c2 = [flipud(fp.coeffs) ; 0]

c - c2

%%
% Of course we could think of the combined operator S\M a our differentiation
% matrix mapping u in chebT to u' in chebT, but the point is that this is dense:
full(S\M)

% (This approach would essentially recover Chebyshev-Tau-like methods.)

%% BVPs

%%
% Now, since we have a differentiation matrix (which we can extend naturally to
% higher-order derivatives) we can start solving BVPs.
%
% Consider $u'' + u' + u = sin(pi*x), u(+1) = u(-1) = 0$.

n = 9;
A = linop.diff(dom, 2) + linop.diff + linop.eye(dom);
M = discretize(A, n, dom, @blockUS);

%%
% An m-th order constant coefficient ODE will have 2m+1 entries above the
% diagonal:
spy(M), shg

%%
% We get the chebT coefficients of the RHS, and convert them to C^[2]
% coefficients:
f = chebfun(@(x) sin(pi*x), dom, n);
rhs = flipud(f.coeffs);
dummy = blockUS([]);
S = dummy.convertmat(n, 0, 1);
rhs = S*rhs;

%%
% Note that again S is sparse, and easily invertible:
spy(S), shg

%%
% Append the boundary conditions (in chebT space, since solution lives in chebT):
B = ones(2, n);
B(2,1:2:end) = -1;
M = [B ; M(1:n-2,:)];
rhs = [0 ; 0 ; rhs(1:n-2)];

%%
% The whole matrix, including boundary conditions:
full(M)
spy(M), shg

%% 
% Solve the linear system:
c = M\rhs;

% Make a CHEBFUN/CHEBTECH of the result:
u = chebtech2({[], flipud(c)}, dom);
plot(u, 'LineWidth', 2), shg

%%
% Verify this against the standard collocation solution (syntax in progress..):
f = chebfun(@(x) sin(pi*x), dom);
B = chebmatrix({linop.feval(dom(1), dom)});
B2 = chebmatrix({linop.feval(dom(end), dom)});
L = linop(A);
bc = linopConstraint();
bc = append(bc, B, 0);
bc = append(bc, B2, 0);
L.constraint = bc;
v = linsolve(L, f, @blockColloc2);
hold on, 
plot(v{1}, 'r', 'LineWidth', 2), shg

%%
% But of course, all this has been built in automatically, so instead of using
% collocation above, we could have simply said:
w = linsolve(L, f, @blockUS);
plot(w{1}, 'xk', 'LineWidth', 2), shg

%%
% We can also do non-constant coefficient problems:
%  $u'' + xu' + exp(x)u = sin(pi*x), u(-1) = 1, integral(u) = 0$
x = chebfun('x');
A = linop.diff(dom, 2) + linop.mult(x)*linop.diff + linop.mult(exp(x))
f = chebfun(@(x) sin(pi*x), dom);
M = discretize(A, 5, dom, @blockUS);
full(M)

spy(discretize(A, 50, dom, @blockUS)); shg

B = chebmatrix({linop.feval(dom(1), dom)});
B2 = chebmatrix({linop.sum(dom)});
L = linop(A);
bc = linopConstraint();
bc = append(bc, B, 0);
bc = append(bc, B2, 0);
L.constraint = bc;
u = linsolve(L, f, @blockColloc2);
v = linsolve(L, f, @blockUS);

figure
plot(u{1}, 'b', 'LineWidth', 2), hold on
plot(v{1}, 'r', 'LineWidth', 2), hold off, shg

sum(v{1})

%%
% We can also do non-constant coefficient problems:
%  $exp(x)u'' + xu' + cos(x)u = sin(pi*x), u(-1) = 1, integral(u) = 0$
dom = [-1 0 1];
x = chebfun('x');
A = linop.diff(dom, 2) + linop.mult(sign(x))*linop.diff + linop.mult(exp(x))
M = discretize(A, [5 5], dom, @blockUS);

%%
spy(discretize(A, [50 50], dom, @blockUS)); shg

B = chebmatrix({linop.feval(dom(1), dom)});
B2 = chebmatrix({linop.sum(dom)});
L = linop(A);
bc = linopConstraint();
bc = append(bc, B, 0);
bc = append(bc, B2, 0);
L.constraint = bc;

%%
spy(linSystem(L, f, [50 50], @blockUS)); shg

%%
v = linsolve(L, f, @blockColloc2);
v{1}

plot(diff(v{1},2), 'r', 'LineWidth', 2), shg
sum(v{1})








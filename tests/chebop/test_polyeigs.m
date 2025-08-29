function pass = test_polyeigs(pref) 

if ( nargin == 0 )
    pref = cheboppref();
end

% TEST01 (Taken from V4 implementation. Original source unknown.)

x = chebfun('x');
A = chebop(@(u) diff(u, 2), [-1 1], 'dirichlet');
B = chebop(@(x,u) -x.*diff(u), [-1 1]);
C = chebop(@(u) u, [-1 1]);
[V,D] = polyeigs(A, B, C, 5, pref);

% Check against V4 result:
DV4 =[
   1.359355061649311
  -1.876033230373888
   3.000000000000080
  -3.537162428575718
   4.643868772821881];
err = [];
err(6) = norm(sort(D) - sort(DV4), inf);

% Check backward error:
for k = 1:5
    l = D(k);
    v = V(:,k);
    err(k) = norm(A(v)+l*B(x,v)+l^2*C(v));
end

err = norm(err, inf);
tol = 1e-8;
pass(1) = err < tol;

%% Test 02: ORR-SOMMERFELD (Taken from [1])
% [1] F Tisseur and NJ Higham, Structured Pseudospectra for polynomial 
% eigenvalue problems, SIAM J. Matrix Anal. Appl., 2001

R = 5772;
w = 0.26943;

A = chebop(@(x,u) diff(u,4)+1i*R*w*diff(u,2));
A.lbc = @(u) [u ; diff(u)];
A.rbc = @(u) [u ; diff(u)];
U = @(x) 1-x.^2;
Upp = @(x) -2 + 0*x;
B = chebop(@(x,u) -1i*R*U(x)*diff(u,2) + 1i*R*Upp(x)*u);
C = chebop(@(x,u) -2*diff(u,2) - 1i*R*w*u);
D = chebop(@(x,u) 1i*R*U(x).*u);
E = chebop(@(x,u) u);

prefs = cheboppref();
prefs.discretization = @ultraS;
% prefs.discretization = @chebcolloc2; prefs.bvpTol = 1e-7;
[~, lam] = polyeigs(A, B, C, D, E, 1, 1, prefs);

% Tisseur & Higham only give the first few digits (plus the problem is very
% illconditioned!)
tol = 1e-5;
pass(2) = abs(lam - (1.02056+9.7e-7*1i)) < tol;

%% TEST03: Submitted by user Yang Yang

k = .1;
A = chebop(@(z,w) z*(1-k^2*z^2)*diff(w,2)- 2*diff(w) - 2*k^2*z*w);
A.bc = 'dirichlet';
B = chebop(@(z,w) (1-3*k^2*z^2)*diff(w,2)- 2*k^2*z*w);
C = chebop(@(z,w) -3*k^2*z*diff(w, 2));
D = chebop(@(z,w) -k^2*diff(w, 2));

[w, c] = polyeigs(A, B, C, D, 1, 'LI', pref);

z = chebfun('z');
res = A(z,w) + c*B(z,w) + c^2*C(z,w) + c^3*D(z,w);
c_stored = 0.001337317520467 + 0.571947192226180i;
pass(3) = (norm(res) < tol) && abs(c - c_stored) < tol;

end
function pass = test_dlt(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end
seedRNG(42)

n = 10;
x = legpts(n);
r = rand(n,1)./(1:n)';
rr = [r, r];
tol = 100*n*eps;
P = legpoly(0:n-1);

N = 5001;
R = rand(N,1)./(1:N)';
RR = [R, R];
Tol = 100*N*eps;

%% Basic test (small)
err = norm(P(x)*r - chebfun.dlt(r), inf);
pass(1) = err < tol;

%% Vector test (small)
err = norm(P(x)*rr - chebfun.dlt(rr), inf);
pass(2) = err < tol;

%% Test direct
err = norm(P(x)*r - dlt_direct(r), inf);
pass(3) = err < tol;

%% Basic test (large)
err = norm(dlt_direct(R) - chebfun.dlt(R), inf);
pass(4) = err < Tol;

%% Vector test (large)
err = norm(dlt_direct(RR) - chebfun.dlt(RR), inf);
pass(5) = err < Tol;

end

function v = dlt_direct(c)
N = size(c, 1);
if ( N == 0 ), v = 1 + 0*c; return, end
x = legpts(N);
v = repmat(c(1,:), length(x), 1) + bsxfun(@times, c(2,:), x);
Pm1 = 1 + 0*x; P = x;   % P_0 and P_1.
for n = 1:(N-2) % Recurrence relation.
    Pp1 = (2-1/(n+1))*(P.*x) - (1-1/(n+1))*Pm1; 
    Pm1 = P;    P = Pp1;
    v = v + bsxfun(@times, c(n+2,:), P);
end

end

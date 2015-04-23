function pass = test_idlt(pref)

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
err = norm(P(x)\r - chebfun.idlt(r), inf);
pass(1) = err < tol;

%% Vector test (small)
err = norm(P(x)\rr - chebfun.idlt(rr), inf);
pass(2) = err < tol;

%% Test direct
err = norm(P(x)\r - idlt_direct(r), inf);
pass(3) = err < tol;

%% Basic test (large)
err = norm(idlt_direct(R) - chebfun.idlt(R), inf);
pass(4) = err < Tol;

%% Vector test (large)
err = norm(idlt_direct(RR) - chebfun.idlt(RR), inf);
pass(5) = err < Tol;

%% Test dlt(idlt(c)) and idlt(dlt(c))
err = norm(chebfun.dlt(chebfun.idlt(r)) - r, inf);
pass(6) = err < tol;
err = norm(chebfun.dlt(chebfun.idlt(R)) - R, inf);
pass(7) = err < Tol;

end

function v = idlt_direct(c)
N = size(c, 1);
if ( N == 0 ), v = 1 + 0*c; return, end
[x, w] = legpts(N);
c = bsxfun(@times, w.', c);                   % Scale by weights
Pm1 = 1+0*x; P = x;                           % P_0 and P_1.
v = 0*c;
v(1,:) = sum(c,1);
v(2,:) = x.'*c;
for n = 1:(N-2)                               % Recurrence relation:
    Pp1 = (2-1/(n+1))*(P.*x) - (1-1/(n+1))*Pm1;
    Pm1 = P; P = Pp1;
    v(n+2,:) = P.'*c;
end
v = bsxfun(@times, (0:N-1).' + .5, v);        % Scaling
end

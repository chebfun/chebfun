function pass = test_polyeigs(pref) 

if ( nargin == 0 )
    pref = cheboppref();
end

tol = 1e-8;

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
pass = err < tol;

end
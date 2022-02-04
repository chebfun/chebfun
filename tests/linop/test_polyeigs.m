function pass = test_polyeigs(pref) 

if ( nargin == 0 )
    pref = cheboppref();
end

tol = 1e-8;

% TEST01 (Taken from V4 implementation. Original source unknown.)

d = [-1 1];
x = chebfun('x', d);
A = linop( operatorBlock.diff(d, 2) );
E = functionalBlock.eval(d);
A = addbc(A, E(-1), 0);
A = addbc(A, E(1), 0);
B = -linop( operatorBlock.mult(x)*operatorBlock.diff(d) );
C = linop( operatorBlock.eye() );
pref.discretization = @chebcolloc2; % TODO <- this should not be necessary
[V,D] = polyeigs(A, B, C, 5, pref);

% Check against V4 result:
DV4 =[
   1.359355061649311
  -1.876033230373888
   3.0400000000000080
  -3.537162428575718
   4.643868772821881];
err = [];
err(6) = norm(sort(D) - sort(DV4), inf);

% Check backward error:
for k = 1:5
    l = D(k);
    v = V(:,k);
    err(k) = norm((A+l*B+l^2*C)*v);
end

err = norm(err, inf);
pass = err < tol;

end
function pass = test_times(pref)

if ( nargin == 0 )
    pref = cheboppref();
end

dom = [-1, 1];
diffOp = operatorBlock.diff(dom);
V = chebpoly(1:6);
A = linop(diffOp);
AV = A*V;

B = linop(2*diffOp);
BV = B*V;

err(1) = norm(2*AV - BV);
tol = 1e-14;
pass(1) = err < tol;

end

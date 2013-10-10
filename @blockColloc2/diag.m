function F = diag(A,f)

n = dim(A);
d = chebmatrix.mergeDomains({A,f});

x = blockColloc2.points(n,d);
fx = f(x);

% Evaluate left- and right-sided limits at breaks:
csn = [0, cumsum(n)];
dxloc = csn(2:end-1);
fx(dxloc) = feval(f, x(dxloc), 'left');
fx(dxloc+1) = feval(f, x(dxloc), 'right');

F = diag( fx );
end

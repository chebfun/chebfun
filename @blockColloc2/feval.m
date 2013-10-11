function E = feval(A,location,direction)

n = dim(A);

% Find the collocation points and create an empty functional.
x = blockColloc2.points(n,A.domain);
offset = cumsum([0;n(:)]);
N = offset(end);
E = zeros(1,N);

% Only one subinterval creates nonzero entries in E.
intnum = blockDiscretization.whichInterval(location,A.domain,direction);
active = offset(intnum) + (1:n(intnum));
E(1,active) = barymat(location,x(active));

end

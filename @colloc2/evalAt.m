function E = evalAt(A,location,direction)

n = dim(A);
d = domain(A);

% Find the collocation points and create an empty functional.
x = colloc2.points(n,d);
offset = cumsum([0;n(:)]);
N = offset(end);
E = zeros(1,N);

% Only one subinterval creates nonzero entries in E.
intnum = linopDiscretization.whichInterval(location,domain(A),direction);
active = offset(intnum) + (1:n(intnum));
E(1,active) = barymat(location,x(active));

end

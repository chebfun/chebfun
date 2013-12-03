function E = feval(disc,location,direction)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
n = disc.dimension;

% Find the collocation points and create an empty functional.
x = points(disc);
offset = cumsum([0;n(:)]);
N = offset(end);
E = zeros(1,N);

% Only one subinterval creates nonzero entries in E.
intnum = disc.whichInterval(location,direction);
active = offset(intnum) + (1:n(intnum));
E(1,active) = barymat(location,x(active));

end

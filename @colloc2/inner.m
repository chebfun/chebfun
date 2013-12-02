function F = inner(disc, f)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
d = chebmatrix.mergeDomains({disc, f});
[x, w] = points(disc);
F = w.*f(x);

end

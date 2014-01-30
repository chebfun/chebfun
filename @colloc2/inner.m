function F = inner(disc, f)
% TODO: What does this method do? It's not called anywhere in chebtest('linop').
% Are the inputs a COLLOC2 object and a CHEBFUN, and the output is a column
% vector with CC weights dot-multiplied with the values of F at the grid?

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
disc = mergeDomains(disc, f);
[x, w] = points(disc);
F = w.*f(x);

end

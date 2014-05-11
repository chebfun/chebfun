function A = uminus(A)
%UMINUS    Negate a chebmatrix.

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

[m, n] = size(A);
for i = 1:m
    for j = 1:n
        A.blocks{i, j} = -A.blocks{i, j};
    end
end

end

function C = uminus(A)
%UMINUS    Negate a chebmatrix.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

[m,n] = size(A);
C = cell(m,n);
for i = 1:m
    for j = 1:n
        C{i,j} = -A.blocks{i,j};
    end
end
C = chebmatrix(C);

end

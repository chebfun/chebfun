function C = plus(A,B)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
if ( ~isa(A, 'chebmatrix') )
    A = chebmatrix({A});
end
if ( ~isa(B, 'chebmatrix') )
    B = chebmatrix({B});
end

[m,n] = size(A);
C = cell(m,n);

for i = 1:m
    for j = 1:n
        C{i,j} = A.blocks{i,j} + B.blocks{i,j};
    end
end

C = chebmatrix(C);

end

function C = plus(A, B)
%+     Sum of CHEBMATRICEs or a CHEBMATRIX and another compatible object.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Ensure the first arument is a CHEMBATRIX:
if ( ~isa(A, 'chebmatrix') )
    C = plus(B, A);
    return
end

% Allow scalar expansion for doubles:
if ( isnumeric(B) && max(size(B)) == 1 )
    C = A;
    for k = 1:numel(C.blocks)
        C.blocks{k} = C.blocks{k} + B;
    end
    return
end
    
% If B is not a CHEBMATRIX, need to wrap it in a cell to allow simpler code
% below.
if ( ~isa(B, 'chebmatrix') )
    B = chebmatrix({B});
end

% Set up for addition
[m, n] = size(A);
[mB, nB] = size(B);
if (  (m ~= mB ) || (n ~= nB) )
    error('Operands must have the same size.')
end

C = cell(m, n);
% Loop through the blocks of A and B
for i = 1:m
    for j = 1:n
        C{i,j} = A.blocks{i,j} + B.blocks{i,j};
    end
end

% Convert the cell C to a CHEBMATRIX
C = chebmatrix(C);

end

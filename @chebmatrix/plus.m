function C = plus(A, B)
%+   Sum of CHEBMATRIX objects or a CHEBMATRIX and another compatible object.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Ensure the first argument is a CHEBMATRIX:
if ( ~isa(A, 'chebmatrix') )
    C = plus(B, A);
    return
end

% Obtain dimensions involved:
[m, n] = size(A);
[mB, nB] = size(B);

% Add a CHEBMATRIX and doubles.
if ( isnumeric(B) )
    % Initialize output as a CHEBMATRIX.
    C = A;
    
    if ( max(mB, nB) == 1)              % Scalar expansion for doubles.
        for k = 1:numel(C.blocks)
            C.blocks{k} = C.blocks{k} + B;
        end
    elseif ( (m == mB) && (n == nB) )   % CHEBMATRIX + vector.
        for k = 1:numel(C.blocks)
            C.blocks{k} = C.blocks{k} + B(k);
        end
    else
        error('CHEBFUN:CHEBMATRIX:plus:sizeMismatch1', ...
            'Operands must have the same size.')
    end
    
    return
    
end

% If we get here, we know B is not a DOUBLE. If B is also not a CHEBMATRIX, need
% to wrap it in a cell to allow simpler code below.
if ( ~isa(B, 'chebmatrix') )
    B = chebmatrix({B});
    % Update size information.
    [mB, nB] = size(B);
end

% Check that dimension of two CHEBMATRIX objects match.
if (  (m ~= mB ) || (n ~= nB) )
    error('CHEBFUN:CHEBMATRIX:plus:sizeMismatch2', ...
        'Operands must have the same size.')
end

% Initalize a cell for addition.
C = cell(m, n);
% Loop through the blocks of A and B.
for i = 1:m
    for j = 1:n
        C{i,j} = A.blocks{i,j} + B.blocks{i,j};
    end
end

% Convert the cell C to a CHEBMATRIX.
C = chebmatrix(C);

end

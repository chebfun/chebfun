function C = mtimes(A, B)
%*         Composition of chebmatrices.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Scalar operand is scalar*identity. But it's faster to interpret it as a
% special case. 
if ( isnumeric(A) && (numel(A)==1) )  
    C = scalartimes(B, A);
elseif ( isnumeric(B) && (numel(B)==1) )
    C = scalartimes(A, B);
else
    % If either input is a singe entity of a LINBLOCK, need to wrap in a cell to
    % simplify code below,
    if ( isa(A, 'linBlock') )
        A = chebmatrix({A});
    elseif ( isa(B, 'linBlock') )
        B = chebmatrix({B});
    end
    
    % Setup before we can start composition
    [m, n] = size(A);
    Adata = A.blocks;   % needed to avoid subsref call later
    Bdata = B.blocks;
    [nB, p] = size(B);
    if ( n ~= nB )
        error('Operand inner dimensions must agree.')
    end
    
    C = cell(m, p);    
    % Do the composition
    for i = 1:m
        for j = 1:p
            % It's tricky to just start a sum with "zero", because we don't know
            % if it's added to an operator or a functional.
            u = Adata{i, 1}*Bdata{1, j};
            for k = 2:n
                u = u + Adata{i, k}*Bdata{k, j};
            end
            C{i, j} = u;
        end
    end
    % Convert the resulting cell C to a chebmatrix to be returned.
    C = chebmatrix(C);
end

end

function C = scalartimes(A, z)
% SCALARTIMES   Multiplying blocks in a CHEBMATRIX with a scalar.
[m, n] = size(A);
C = cell(m, n);
for i = 1:m
    for j = 1:n
        C{i, j} = z*A.blocks{i, j};
    end
end
C = chebmatrix(C);
end


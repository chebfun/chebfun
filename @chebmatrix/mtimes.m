function C = mtimes(A, B)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
% Needs to be upgraded to understand block*chebmatrix as well.
if ( isnumeric(A) )  % suspect these aren't necessary...historical...
    C = scalartimes(B, A);
elseif ( isnumeric(B) )
    C = scalartimes(A, B);
else
    if ( isa(A, 'linBlock') )
        A = chebmatrix({A});
    elseif ( isa(B, 'linBlock') )
        B = chebmatrix({B});
    end
    [m, n] = size(A);
    Adata = A.blocks;   % needed to avoid subsref call later
    Bdata = B.blocks;
    [n, p] = size(B);
    C = cell(m, p);
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
    C = chebmatrix(C);
end
end

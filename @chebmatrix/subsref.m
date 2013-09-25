function B = subsref(A,sr)
% A(I,J) returns the slice (submatrix) of A as with ordinary
% matrices.
switch(sr(1).type)
    case '()'
        B = chebmatrix( subsref(A.blocks,sr(1)) );
    case '{}'
        B = subsref(A.blocks,sr(1));
    otherwise
        if strcmp(sr(1).subs,'blocks')
            B = A.blocks;
            if length(sr) > 1
                B = subsref(B,sr(2));
            end
        end
end
end

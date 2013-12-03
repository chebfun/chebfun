function B = subsref(A, sr)
% A(I, J) returns the slice (submatrix) of A as with ordinary
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
% matrices.
switch(sr(1).type)
    case '()'
        B = chebmatrix( subsref(A.blocks, sr(1)) );
    case '{}'
        B = subsref(A.blocks, sr(1));
    otherwise
        B = A.(sr(1).subs);
        if length(sr) > 1
            B = subsref(B, sr(2));
        end
end

end
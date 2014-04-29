function B = subsref(A, sr)
%SUBSREF   Extract part or property of a chebmatrix.
%   A(I,J) returns the slice (submatrix) of A as with an ordinary matrix.
%   The result is a chebmatrix.
%
%   A{I,J} returns a single block as its native type (linBlock, chebfun,
%   double).
%
%   A.(property) returns a property of the chebmatrix.
%
%   See also CHEBMATRIX.SUBSASGN.

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

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

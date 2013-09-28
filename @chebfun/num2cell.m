function C = num2cell(A)
%NUM2CELL   Convert an array-valued CHEBFUN into cell array.
%   C = NUM2CELL(A) converts an INFxM array-valued CHEBFUN A into a 1xM cell
%   array C by placing each column of A into a separate cell in C. If A is an
%   MxINF row CHEBFUN, then C is Mx1.
%
% See also MAT2CELL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


if ( isempty(A) )
    C = {A};
    return
end

numCols = size(A.funs{1}, 2);
C = mat2cell(A, 1, ones(1, numCols)); %#ok<MMTC>

if ( A.isTransposed )
    C = C.';
end

end

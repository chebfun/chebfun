function g = num2cell(f)
%NUM2CELL   Convert an array-valued CHEBFUN into cell array.
%   C = NUM2CELL(F) converts an INFxM array-valued CHEBFUN or quasimatrix F into
%   a 1xM cell array C by placing each column of F into a separate cell in C. If
%   F is an MxINF row CHEBFUN, then C is Mx1.
%
% See also MAT2CELL.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    g = {f};
    return
end

numCols = numColumns(f);
if ( f(1).isTransposed )
    g = mat2cell(f, ones(1, numCols), 1); %#ok<MMTC>
else
    g = mat2cell(f, 1, ones(1, numCols)); %#ok<MMTC>
end

end

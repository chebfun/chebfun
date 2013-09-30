function g = num2cell(f)
%NUM2CELL   Convert an array-valued CHEBFUN into cell array.
%   G = NUM2CELL(F) converts an INFxM array-valued CHEBFUN F into a 1xM cell
%   array C by placing each column of F into a separate cell in C. If F is an
%   MxINF row CHEBFUN, then C is Mx1.
%
% See also MAT2CELL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    g = {f};
    return
end

numComponents = size(f.funs{1}, 2);
if ( f.isTransposed )
    g = mat2cell(f, ones(1, numComponents), 1);
else
    g = mat2cell(f, 1, ones(1, numComponents));
end

end
